/*  Lzlib - A compression library for lzip files
    Copyright (C) 2009, 2010 Antonio Diaz Diaz.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.

    As a special exception, you may use this file as part of a free
    software library without restriction.  Specifically, if other files
    instantiate templates or use macros or inline functions from this
    file, or you compile this file and link it with other files to
    produce an executable, this file does not by itself cause the
    resulting executable to be covered by the GNU General Public
    License.  This exception does not however invalidate any other
    reasons why the executable file might be covered by the GNU General
    Public License.
*/

class Input_buffer : public Circular_buffer
  {
  enum { min_available_bytes = 8 };
  bool at_stream_end_;

public:
  Input_buffer()
    :
    Circular_buffer( 65536 + min_available_bytes ),
    at_stream_end_( false ) {}

  bool at_stream_end() const throw() { return at_stream_end_; }
  void finish() throw() { at_stream_end_ = true; }
  bool finished() const throw() { return at_stream_end_ && !used_bytes(); }
  int free_bytes() const throw()
    { if( at_stream_end_ ) return 0; return Circular_buffer::free_bytes(); }
  void purge() throw() { at_stream_end_ = true; Circular_buffer::reset(); }
  void reset() throw() { at_stream_end_ = false; Circular_buffer::reset(); }

  bool find_header() throw();
  bool read_header( File_header & header ) throw();

  bool enough_available_bytes() const throw()
    {
    return ( used_bytes() > 0 &&
           ( at_stream_end_ || used_bytes() >= min_available_bytes ) );
    }

  int write_data( const uint8_t * const in_buffer, const int in_size ) throw()
    {
    if( at_stream_end_ || in_size <= 0 ) return 0;
    return Circular_buffer::write_data( in_buffer, in_size );
    }
  };


class Range_decoder
  {
  mutable long long member_pos;
  uint32_t code;
  uint32_t range;
  bool reload_pending;
  Input_buffer & ibuf;

public:
  Range_decoder( const int header_size, Input_buffer & buf )
    :
    member_pos( header_size ),
    code( 0 ),
    range( 0xFFFFFFFF ),
    reload_pending( false ),
    ibuf( buf )
    { for( int i = 0; i < 5; ++i ) code = (code << 8) | get_byte(); }

  bool at_stream_end() const throw() { return ibuf.at_stream_end(); }
  int available_bytes() const throw() { return ibuf.used_bytes(); }
  bool code_is_zero() const throw() { return ( code == 0 ); }
  bool enough_available_bytes() const throw()
    { return ibuf.enough_available_bytes(); }
  bool finished() const throw() { return ibuf.finished(); }
  long long member_position() const throw() { return member_pos; }

  uint8_t get_byte() const
    {
    ++member_pos;
    return ibuf.get_byte();
    }

  bool try_reload( const bool force = false ) throw()
    {
    if( force ) reload_pending = true;
    if( reload_pending && available_bytes() >= 5 )
      {
      reload_pending = false;
      code = 0;
      range = 0xFFFFFFFF;
      for( int i = 0; i < 5; ++i ) code = (code << 8) | get_byte();
      }
    return !reload_pending;
    }

  void normalize()
    {
    if( range <= 0x00FFFFFF )
      { range <<= 8; code = (code << 8) | get_byte(); }
    }

  int decode( const int num_bits )
    {
    int symbol = 0;
    for( int i = num_bits; i > 0; --i )
      {
      symbol <<= 1;
      if( range <= 0x00FFFFFF )
        {
        range <<= 7; code = (code << 8) | get_byte();
        if( code >= range ) { code -= range; symbol |= 1; }
        }
      else
        {
        range >>= 1;
        if( code >= range ) { code -= range; symbol |= 1; }
        }
      }
    return symbol;
    }

  int decode_bit( Bit_model & bm )
    {
    normalize();
    const uint32_t bound = ( range >> bit_model_total_bits ) * bm.probability;
    if( code < bound )
      {
      range = bound;
      bm.probability += (bit_model_total - bm.probability) >> bit_model_move_bits;
      return 0;
      }
    else
      {
      range -= bound;
      code -= bound;
      bm.probability -= bm.probability >> bit_model_move_bits;
      return 1;
      }
    }

  int decode_tree( Bit_model bm[], const int num_bits )
    {
    int model = 1;
    for( int i = num_bits; i > 0; --i )
      model = ( model << 1 ) | decode_bit( bm[model] );
    return model - (1 << num_bits);
    }

  int decode_tree_reversed( Bit_model bm[], const int num_bits )
    {
    int model = 1;
    int symbol = 0;
    for( int i = 0; i < num_bits; ++i )
      {
      const int bit = decode_bit( bm[model] );
      model <<= 1;
      if( bit ) { model |= 1; symbol |= (1 << i); }
      }
    return symbol;
    }

  int decode_matched( Bit_model bm[], const int match_byte )
    {
    Bit_model *bm1 = bm + 0x100;
    int symbol = 1;
    for( int i = 1; i <= 8; ++i )
      {
      const int match_bit = ( match_byte << i ) & 0x100;
      const int bit = decode_bit( bm1[match_bit+symbol] );
      symbol = ( symbol << 1 ) | bit;
      if( ( match_bit && !bit ) || ( !match_bit && bit ) )
        {
        while( ++i <= 8 )
          symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
        break;
        }
      }
    return symbol & 0xFF;
    }
  };


class Len_decoder
  {
  Bit_model choice1;
  Bit_model choice2;
  Bit_model bm_low[pos_states][len_low_symbols];
  Bit_model bm_mid[pos_states][len_mid_symbols];
  Bit_model bm_high[len_high_symbols];

public:
  int decode( Range_decoder & range_decoder, const int pos_state )
    {
    if( range_decoder.decode_bit( choice1 ) == 0 )
      return range_decoder.decode_tree( bm_low[pos_state], len_low_bits );
    if( range_decoder.decode_bit( choice2 ) == 0 )
      return len_low_symbols +
             range_decoder.decode_tree( bm_mid[pos_state], len_mid_bits );
    return len_low_symbols + len_mid_symbols +
           range_decoder.decode_tree( bm_high, len_high_bits );
    }
  };


class Literal_decoder
  {
  Bit_model bm_literal[1<<literal_context_bits][0x300];

  int state( const int prev_byte ) const throw()
    { return ( prev_byte >> ( 8 - literal_context_bits ) ); }

public:
  uint8_t decode( Range_decoder & range_decoder, const uint8_t prev_byte )
    { return range_decoder.decode_tree( bm_literal[state(prev_byte)], 8 ); }

  uint8_t decode_matched( Range_decoder & range_decoder,
                          const uint8_t prev_byte, const uint8_t match_byte )
    { return range_decoder.decode_matched( bm_literal[state(prev_byte)], match_byte ); }
  };


class LZ_decoder : public Circular_buffer
  {
  enum { min_free_bytes = max_match_len };
  long long partial_data_pos;
  const int format_version;
  const int dictionary_size;
  uint32_t crc_;
  bool member_finished_;
  bool verify_trailer_pending;
  unsigned int rep0;		// rep[0-3] latest four distances
  unsigned int rep1;		// used for efficient coding of
  unsigned int rep2;		// repeated distances
  unsigned int rep3;
  State state;

  Bit_model bm_match[State::states][pos_states];
  Bit_model bm_rep[State::states];
  Bit_model bm_rep0[State::states];
  Bit_model bm_rep1[State::states];
  Bit_model bm_rep2[State::states];
  Bit_model bm_len[State::states][pos_states];
  Bit_model bm_dis_slot[max_dis_states][1<<dis_slot_bits];
  Bit_model bm_dis[modeled_distances-end_dis_model];
  Bit_model bm_align[dis_align_size];

  Range_decoder range_decoder;
  Len_decoder len_decoder;
  Len_decoder rep_match_len_decoder;
  Literal_decoder literal_decoder;

  uint8_t get_byte( const int distance ) const throw()
    {
    int i = put - distance - 1;
    if( i < 0 ) i += buffer_size;
    return buffer[i];
    }

  void put_byte( const uint8_t b )
    {
    crc32.update( crc_, b );
    buffer[put] = b;
    if( ++put >= buffer_size ) { partial_data_pos += put; put = 0; }
    }

  void copy_block( const int distance, int len )
    {
    int i = put - distance - 1;
    if( i < 0 ) i += buffer_size;
    if( len < buffer_size - std::max( put, i ) && len <= std::abs( put - i ) )
      {
      crc32.update( crc_, buffer + i, len );
      std::memcpy( buffer + put, buffer + i, len );
      put += len;
      }
    else for( ; len > 0 ; --len )
      {
      crc32.update( crc_, buffer[i] );
      buffer[put] = buffer[i];
      if( ++put >= buffer_size ) { partial_data_pos += put; put = 0; }
      if( ++i >= buffer_size ) i = 0;
      }
    }

  bool verify_trailer();

public:
  LZ_decoder( const File_header & header, Input_buffer & ibuf )
    :
    Circular_buffer( std::max( 65536, header.dictionary_size() ) + min_free_bytes ),
    partial_data_pos( 0 ),
    format_version( header.version ),
    dictionary_size( header.dictionary_size() ),
    crc_( 0xFFFFFFFF ),
    member_finished_( false ),
    verify_trailer_pending( false ),
    rep0( 0 ),
    rep1( 0 ),
    rep2( 0 ),
    rep3( 0 ),
    range_decoder( sizeof header, ibuf ),
    literal_decoder()
    { buffer[buffer_size-1] = 0; }	// prev_byte of first_byte

  bool enough_free_bytes() const throw()
    { return free_bytes() >= min_free_bytes; }

  uint32_t crc() const throw() { return crc_ ^ 0xFFFFFFFF; }
  int decode_member();
  bool member_finished() const throw()
    { return ( member_finished_ && !used_bytes() ); }

  long long member_position() const throw()
    { return range_decoder.member_position(); }
  long long data_position() const throw()
    { return partial_data_pos + put; }
  };
