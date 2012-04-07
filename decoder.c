/*  Lzlib - A compression library for lzip files
    Copyright (C) 2009, 2010, 2011, 2012 Antonio Diaz Diaz.

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

struct Circular_buffer
  {
  uint8_t * buffer;
  int buffer_size;		/* capacity == buffer_size - 1 */
  int get;			/* buffer is empty when get == put */
  int put;
  };

static inline void Cb_reset( struct Circular_buffer * const cb )
  { cb->get = 0; cb->put = 0; }

static inline bool Cb_init( struct Circular_buffer * const cb,
                            const int buf_size )
  {
  cb->buffer = (uint8_t *)malloc( buf_size + 1 );
  cb->buffer_size = buf_size + 1;
  cb->get = 0;
  cb->put = 0;
  return ( cb->buffer != 0 );
  }

static inline void Cb_free( struct Circular_buffer * const cb )
  { free( cb->buffer ); cb->buffer = 0; }

static inline int Cb_used_bytes( const struct Circular_buffer * const cb )
    { return ( (cb->get <= cb->put) ? 0 : cb->buffer_size ) + cb->put - cb->get; }

static inline int Cb_free_bytes( const struct Circular_buffer * const cb )
    { return ( (cb->get <= cb->put) ? cb->buffer_size : 0 ) - cb->put + cb->get - 1; }

static inline uint8_t Cb_get_byte( struct Circular_buffer * const cb )
    {
    const uint8_t b = cb->buffer[cb->get];
    if( ++cb->get >= cb->buffer_size ) cb->get = 0;
    return b;
    }

static inline void Cb_put_byte( struct Circular_buffer * const cb,
                                const uint8_t b )
    {
    cb->buffer[cb->put] = b;
    if( ++cb->put >= cb->buffer_size ) cb->put = 0;
    }


/* Copies up to 'out_size' bytes to 'out_buffer' and updates 'get'.
   Returns the number of bytes copied.
*/
static int Cb_read_data( struct Circular_buffer * const cb,
                         uint8_t * const out_buffer, const int out_size )
  {
  if( out_size < 0 ) return 0;
  int size = 0;
  if( cb->get > cb->put )
    {
    size = min( cb->buffer_size - cb->get, out_size );
    if( size > 0 )
      {
      memcpy( out_buffer, cb->buffer + cb->get, size );
      cb->get += size;
      if( cb->get >= cb->buffer_size ) cb->get = 0;
      }
    }
  if( cb->get < cb->put )
    {
    const int size2 = min( cb->put - cb->get, out_size - size );
    if( size2 > 0 )
      {
      memcpy( out_buffer + size, cb->buffer + cb->get, size2 );
      cb->get += size2;
      size += size2;
      }
    }
  return size;
  }


/* Copies up to 'in_size' bytes from 'in_buffer' and updates 'put'.
   Returns the number of bytes copied.
*/
static int Cb_write_data( struct Circular_buffer * const cb,
                          const uint8_t * const in_buffer, const int in_size )
  {
  if( in_size < 0 ) return 0;
  int size = 0;
  if( cb->put >= cb->get )
    {
    size = min( cb->buffer_size - cb->put - (cb->get == 0), in_size );
    if( size > 0 )
      {
      memcpy( cb->buffer + cb->put, in_buffer, size );
      cb->put += size;
      if( cb->put >= cb->buffer_size ) cb->put = 0;
      }
    }
  if( cb->put < cb->get )
    {
    const int size2 = min( cb->get - cb->put - 1, in_size - size );
    if( size2 > 0 )
      {
      memcpy( cb->buffer + cb->put, in_buffer + size, size2 );
      cb->put += size2;
      size += size2;
      }
    }
  return size;
  }


enum { rd_min_available_bytes = 8 };

struct Range_decoder
  {
  struct Circular_buffer cb;
  long long member_position;
  uint32_t code;
  uint32_t range;
  bool reload_pending;
  bool at_stream_end;
  };

static inline bool Rd_init( struct Range_decoder * const rdec )
  {
  if( !Cb_init( &rdec->cb, 65536 + rd_min_available_bytes ) ) return false;
  rdec->member_position = 0;
  rdec->code = 0;
  rdec->range = 0xFFFFFFFFU;
  rdec->reload_pending = false;
  rdec->at_stream_end = false;
  return true;
  }

static inline void Rd_free( struct Range_decoder * const rdec )
  { Cb_free( &rdec->cb ); }

static inline int Rd_available_bytes( const struct Range_decoder * const rdec )
  { return Cb_used_bytes( &rdec->cb ); }

static inline void Rd_finish( struct Range_decoder * const rdec )
  { rdec->at_stream_end = true; }

static inline bool Rd_finished( const struct Range_decoder * const rdec )
  { return rdec->at_stream_end && !Cb_used_bytes( &rdec->cb ); }

static inline int Rd_free_bytes( const struct Range_decoder * const rdec )
  { if( rdec->at_stream_end ) return 0; return Cb_free_bytes( &rdec->cb ); }

static inline void Rd_purge( struct Range_decoder * const rdec )
  { rdec->at_stream_end = true; Cb_reset( &rdec->cb ); }

static inline void Rd_reset( struct Range_decoder * const rdec )
  { rdec->at_stream_end = false; Cb_reset( &rdec->cb ); }


/* Seeks a member header and updates 'get'.
   Returns true if it finds a valid header.
*/
static bool Rd_find_header( struct Range_decoder * const rdec )
  {
  while( rdec->cb.get != rdec->cb.put )
    {
    if( rdec->cb.buffer[rdec->cb.get] == magic_string[0] )
      {
      int g = rdec->cb.get;
      int i;
      File_header header;
      for( i = 0; i < Fh_size; ++i )
        {
        if( g == rdec->cb.put ) return false;	/* not enough data */
        header[i] = rdec->cb.buffer[g];
        if( ++g >= rdec->cb.buffer_size ) g = 0;
        }
      if( Fh_verify( header ) ) return true;
      }
    if( ++rdec->cb.get >= rdec->cb.buffer_size ) rdec->cb.get = 0;
    }
  return false;
  }


/* Returns true, fills 'header', and updates 'get' if 'get' points to a
   valid header.
   Else returns false and leaves 'get' unmodified.
*/
static bool Rd_read_header( struct Range_decoder * const rdec,
                            File_header header )
  {
  int g = rdec->cb.get;
  int i;
  for( i = 0; i < Fh_size; ++i )
    {
    if( g == rdec->cb.put ) return false;	/* not enough data */
    header[i] = rdec->cb.buffer[g];
    if( ++g >= rdec->cb.buffer_size ) g = 0;
    }
  if( Fh_verify( header ) )
    {
    rdec->cb.get = g;
    rdec->member_position = Fh_size;
    rdec->reload_pending = true;
    return true;
    }
  return false;
  }


static inline int Rd_write_data( struct Range_decoder * const rdec,
                                 const uint8_t * const inbuf, const int size )
  {
  if( rdec->at_stream_end || size <= 0 ) return 0;
  return Cb_write_data( &rdec->cb, inbuf, size );
  }

static inline uint8_t Rd_get_byte( struct Range_decoder * const rdec )
  {
  ++rdec->member_position;
  return Cb_get_byte( &rdec->cb );
  }


static bool Rd_try_reload( struct Range_decoder * const rdec, const bool force )
  {
  if( force ) rdec->reload_pending = true;
  if( rdec->reload_pending && Rd_available_bytes( rdec ) >= 5 )
    {
    int i;
    rdec->reload_pending = false;
    rdec->code = 0;
    rdec->range = 0xFFFFFFFFU;
    for( i = 0; i < 5; ++i )
      rdec->code = (rdec->code << 8) | Rd_get_byte( rdec );
    }
  return !rdec->reload_pending;
  }


static inline int Rd_decode( struct Range_decoder * const rdec,
                             const int num_bits )
  {
  int symbol = 0;
  int i;
  for( i = num_bits; i > 0; --i )
    {
    symbol <<= 1;
    if( rdec->range <= 0x00FFFFFFU )
      {
      rdec->range <<= 7;
      rdec->code = (rdec->code << 8) | Rd_get_byte( rdec );
      if( rdec->code >= rdec->range )
        { rdec->code -= rdec->range; symbol |= 1; }
      }
    else
      {
      rdec->range >>= 1;
      if( rdec->code >= rdec->range )
        { rdec->code -= rdec->range; symbol |= 1; }
      }
    }
  return symbol;
  }


static inline void Rd_normalize( struct Range_decoder * const rdec )
  {
  if( rdec->range <= 0x00FFFFFFU )
    {
    rdec->range <<= 8;
    rdec->code = (rdec->code << 8) | Rd_get_byte( rdec );
    }
  }


static inline int Rd_decode_bit( struct Range_decoder * const rdec,
                                 Bit_model * const probability )
  {
  uint32_t bound;
  Rd_normalize( rdec );
  bound = ( rdec->range >> bit_model_total_bits ) * *probability;
  if( rdec->code < bound )
    {
    rdec->range = bound;
    *probability += (bit_model_total - *probability) >> bit_model_move_bits;
    return 0;
    }
  else
    {
    rdec->range -= bound;
    rdec->code -= bound;
    *probability -= *probability >> bit_model_move_bits;
    return 1;
    }
  }


static inline int Rd_decode_tree( struct Range_decoder * const rdec,
                                  Bit_model bm[], const int num_bits )
  {
  int model = 1;
  int i;
  for( i = num_bits; i > 0; --i )
    model = ( model << 1 ) | Rd_decode_bit( rdec, &bm[model] );
  return model - (1 << num_bits);
  }


static inline int Rd_decode_matched( struct Range_decoder * const rdec,
                                     Bit_model bm[], const int match_byte )
  {
  Bit_model * const bm1 = bm + 0x100;
  int symbol = 1;
  int i;
  for( i = 7; i >= 0; --i )
    {
    const int match_bit = ( match_byte >> i ) & 1;
    const int bit = Rd_decode_bit( rdec, &bm1[(match_bit<<8)+symbol] );
    symbol = ( symbol << 1 ) | bit;
    if( match_bit != bit )
      {
      while( --i >= 0 )
        symbol = ( symbol << 1 ) | Rd_decode_bit( rdec, &bm[symbol] );
      break;
      }
    }
  return symbol & 0xFF;
  }


static inline int Rd_decode_tree_reversed( struct Range_decoder * const rdec,
                                           Bit_model bm[], const int num_bits )
  {
  int model = 1;
  int symbol = 0;
  int i;
  for( i = 0; i < num_bits; ++i )
    {
    const int bit = Rd_decode_bit( rdec, &bm[model] );
    model <<= 1;
    if( bit ) { model |= 1; symbol |= (1 << i); }
    }
  return symbol;
  }


static inline bool Rd_enough_available_bytes( const struct Range_decoder * const rdec )
  {
  return ( Cb_used_bytes( &rdec->cb ) >= rd_min_available_bytes ||
           ( rdec->at_stream_end && Cb_used_bytes( &rdec->cb ) > 0 ) );
  }


static inline int Rd_read_data( struct Range_decoder * const rdec,
                                uint8_t * const outbuf, const int size )
  {
  const int sz = Cb_read_data( &rdec->cb, outbuf, size );
  if( sz > 0 ) rdec->member_position += sz;
  return sz;
  }


struct Len_decoder
  {
  Bit_model choice1;
  Bit_model choice2;
  Bit_model bm_low[pos_states][len_low_symbols];
  Bit_model bm_mid[pos_states][len_mid_symbols];
  Bit_model bm_high[len_high_symbols];
  };

static inline void Led_init( struct Len_decoder * const len_decoder )
  {
  int i, j;
  Bm_init( &len_decoder->choice1 );
  Bm_init( &len_decoder->choice2 );
  for( i = 0; i < pos_states; ++i )
    for( j = 0; j < len_low_symbols; ++j )
      Bm_init( &len_decoder->bm_low[i][j] );
  for( i = 0; i < pos_states; ++i )
    for( j = 0; j < len_mid_symbols; ++j )
      Bm_init( &len_decoder->bm_mid[i][j] );
  for( i = 0; i < len_high_symbols; ++i )
    Bm_init( &len_decoder->bm_high[i] );
  }

static inline int Led_decode( struct Len_decoder * const len_decoder,
                              struct Range_decoder * const rdec,
                              const int pos_state )
  {
  if( Rd_decode_bit( rdec, &len_decoder->choice1 ) == 0 )
    return Rd_decode_tree( rdec, len_decoder->bm_low[pos_state], len_low_bits );
  if( Rd_decode_bit( rdec, &len_decoder->choice2 ) == 0 )
    return len_low_symbols +
           Rd_decode_tree( rdec, len_decoder->bm_mid[pos_state], len_mid_bits );
  return len_low_symbols + len_mid_symbols +
         Rd_decode_tree( rdec, len_decoder->bm_high, len_high_bits );
  }


struct Literal_decoder
  {
  Bit_model bm_literal[1<<literal_context_bits][0x300];
  };

static inline void Lid_init( struct Literal_decoder * const lidec )
  {
  int i, j;
  for( i = 0; i < 1<<literal_context_bits; ++i )
    for( j = 0; j < 0x300; ++j )
      Bm_init( &lidec->bm_literal[i][j] );
  }

static inline int Lid_state( const uint8_t prev_byte )
  { return ( prev_byte >> ( 8 - literal_context_bits ) ); }

static inline uint8_t Lid_decode( struct Literal_decoder * const lidec,
                                  struct Range_decoder * const rdec,
                                  const uint8_t prev_byte )
  { return Rd_decode_tree( rdec, lidec->bm_literal[Lid_state(prev_byte)], 8 ); }

static inline uint8_t Lid_decode_matched( struct Literal_decoder * const lidec,
                                          struct Range_decoder * const rdec,
                                          const uint8_t prev_byte,
                                          const uint8_t match_byte )
  { return Rd_decode_matched( rdec, lidec->bm_literal[Lid_state(prev_byte)], match_byte ); }


enum { lzd_min_free_bytes = max_match_len };

struct LZ_decoder
  {
  struct Circular_buffer cb;
  long long partial_data_pos;
  int dictionary_size;
  uint32_t crc;
  int member_version;
  bool member_finished;
  bool verify_trailer_pending;
  unsigned int rep0;		/* rep[0-3] latest four distances */
  unsigned int rep1;		/* used for efficient coding of */
  unsigned int rep2;		/* repeated distances */
  unsigned int rep3;
  State state;

  Bit_model bm_match[states][pos_states];
  Bit_model bm_rep[states];
  Bit_model bm_rep0[states];
  Bit_model bm_rep1[states];
  Bit_model bm_rep2[states];
  Bit_model bm_len[states][pos_states];
  Bit_model bm_dis_slot[max_dis_states][1<<dis_slot_bits];
  Bit_model bm_dis[modeled_distances-end_dis_model+1];
  Bit_model bm_align[dis_align_size];

  struct Range_decoder * range_decoder;
  struct Len_decoder len_decoder;
  struct Len_decoder rep_match_len_decoder;
  struct Literal_decoder literal_decoder;
  };

static inline bool LZd_init( struct LZ_decoder * const decoder,
                             const File_header header,
                             struct Range_decoder * const rdec )
  {
  int i, j;
  if( !Cb_init( &decoder->cb, max( 65536, Fh_get_dictionary_size( header ) ) + lzd_min_free_bytes ) )
    return false;
  decoder->partial_data_pos = 0;
  decoder->dictionary_size = Fh_get_dictionary_size( header );
  decoder->crc = 0xFFFFFFFFU;
  decoder->member_version = Fh_version( header );
  decoder->member_finished = false;
  decoder->verify_trailer_pending = false;
  decoder->rep0 = 0;
  decoder->rep1 = 0;
  decoder->rep2 = 0;
  decoder->rep3 = 0;
  decoder->state = 0;

  for( i = 0; i < states; ++i )
    {
    for( j = 0; j < pos_states; ++j )
      {
      Bm_init( &decoder->bm_match[i][j] );
      Bm_init( &decoder->bm_len[i][j] );
      }
    Bm_init( &decoder->bm_rep[i] );
    Bm_init( &decoder->bm_rep0[i] );
    Bm_init( &decoder->bm_rep1[i] );
    Bm_init( &decoder->bm_rep2[i] );
    }
  for( i = 0; i < max_dis_states; ++i )
    for( j = 0; j < 1<<dis_slot_bits; ++j )
      Bm_init( &decoder->bm_dis_slot[i][j] );
  for( i = 0; i < modeled_distances-end_dis_model+1; ++i )
    Bm_init( &decoder->bm_dis[i] );
  for( i = 0; i < dis_align_size; ++i )
    Bm_init( &decoder->bm_align[i] );

  decoder->range_decoder = rdec;
  Led_init( &decoder->len_decoder );
  Led_init( &decoder->rep_match_len_decoder );
  Lid_init( &decoder->literal_decoder );
  decoder->cb.buffer[decoder->cb.buffer_size-1] = 0; /* prev_byte of first_byte */
  return true;
  }

static inline void LZd_free( struct LZ_decoder * const decoder )
  { Cb_free( &decoder->cb ); }

static inline bool LZd_member_finished( const struct LZ_decoder * const decoder )
  { return ( decoder->member_finished && !Cb_used_bytes( &decoder->cb ) ); }

static inline uint32_t LZd_crc( const struct LZ_decoder * const decoder )
  { return decoder->crc ^ 0xFFFFFFFFU; }

static inline long long LZd_data_position( const struct LZ_decoder * const decoder )
  { return decoder->partial_data_pos + decoder->cb.put; }


static bool LZd_verify_trailer( struct LZ_decoder * const decoder )
  {
  File_trailer trailer;
  const int trailer_size = Ft_versioned_size( decoder->member_version );
  const long long member_size =
    decoder->range_decoder->member_position + trailer_size;

  int size = Rd_read_data( decoder->range_decoder, trailer, trailer_size );
  if( size < trailer_size ) return false;

  if( decoder->member_version == 0 ) Ft_set_member_size( trailer, member_size );

  if( decoder->range_decoder->code != 0 ||
      Ft_get_data_crc( trailer ) != LZd_crc( decoder ) ||
      Ft_get_data_size( trailer ) != LZd_data_position( decoder ) ||
      Ft_get_member_size( trailer ) != member_size ) return false;
  return true;
  }


static inline void LZd_copy_block( struct LZ_decoder * const decoder,
                                   const int distance, int len )
  {
  int i = decoder->cb.put - distance - 1;
  if( i < 0 ) i += decoder->cb.buffer_size;
  if( len < decoder->cb.buffer_size - max( decoder->cb.put, i ) &&
      len <= abs( decoder->cb.put - i ) )
    {
    CRC32_update_buf( &decoder->crc, decoder->cb.buffer + i, len );
    memcpy( decoder->cb.buffer + decoder->cb.put, decoder->cb.buffer + i, len );
    decoder->cb.put += len;
    }
  else for( ; len > 0; --len )
    {
    CRC32_update_byte( &decoder->crc, decoder->cb.buffer[i] );
    decoder->cb.buffer[decoder->cb.put] = decoder->cb.buffer[i];
    if( ++decoder->cb.put >= decoder->cb.buffer_size )
      { decoder->partial_data_pos += decoder->cb.put; decoder->cb.put = 0; }
    if( ++i >= decoder->cb.buffer_size ) i = 0;
    }
  }


static inline bool LZd_enough_free_bytes( const struct LZ_decoder * const decoder )
  { return Cb_free_bytes( &decoder->cb ) >= lzd_min_free_bytes; }


static inline uint8_t LZd_get_byte( const struct LZ_decoder * const decoder,
                                    const int distance )
  {
  int i = decoder->cb.put - distance - 1;
  if( i < 0 ) i += decoder->cb.buffer_size;
  return decoder->cb.buffer[i];
  }

static inline uint8_t LZd_get_prev_byte( const struct LZ_decoder * const decoder )
  {
  const int i =
    ( ( decoder->cb.put > 0 ) ? decoder->cb.put : decoder->cb.buffer_size ) - 1;
  return decoder->cb.buffer[i];
  }

static inline void LZd_put_byte( struct LZ_decoder * const decoder,
                                 const uint8_t b )
  {
  CRC32_update_byte( &decoder->crc, b );
  decoder->cb.buffer[decoder->cb.put] = b;
  if( ++decoder->cb.put >= decoder->cb.buffer_size )
    { decoder->partial_data_pos += decoder->cb.put; decoder->cb.put = 0; }
  }


/* Return value: 0 = OK, 1 = decoder error, 2 = unexpected EOF,
                 3 = trailer error, 4 = unknown marker found. */
static int LZd_decode_member( struct LZ_decoder * const decoder )
  {
  State * const state = &decoder->state;
  if( decoder->member_finished ) return 0;
  if( !Rd_try_reload( decoder->range_decoder, false ) ) return 0;
  if( decoder->verify_trailer_pending )
    {
    if( Rd_available_bytes( decoder->range_decoder ) < Ft_versioned_size( decoder->member_version ) &&
        !decoder->range_decoder->at_stream_end )
      return 0;
    decoder->verify_trailer_pending = false;
    decoder->member_finished = true;
    if( LZd_verify_trailer( decoder ) ) return 0; else return 3;
    }

  while( !Rd_finished( decoder->range_decoder ) )
    {
    const int pos_state = LZd_data_position( decoder ) & pos_state_mask;
    if( !Rd_enough_available_bytes( decoder->range_decoder ) ||
        !LZd_enough_free_bytes( decoder ) )
      return 0;
    if( Rd_decode_bit( decoder->range_decoder, &decoder->bm_match[*state][pos_state] ) == 0 )
      {
      const uint8_t prev_byte = LZd_get_prev_byte( decoder );
      if( St_is_char( *state ) )
        LZd_put_byte( decoder, Lid_decode( &decoder->literal_decoder,
          decoder->range_decoder, prev_byte ) );
      else
        LZd_put_byte( decoder, Lid_decode_matched( &decoder->literal_decoder,
          decoder->range_decoder, prev_byte, LZd_get_byte( decoder, decoder->rep0 ) ) );
      St_set_char( state );
      }
    else
      {
      int len;
      if( Rd_decode_bit( decoder->range_decoder, &decoder->bm_rep[*state] ) == 1 )
        {
        len = 0;
        if( Rd_decode_bit( decoder->range_decoder, &decoder->bm_rep0[*state] ) == 1 )
          {
          unsigned int distance;
          if( Rd_decode_bit( decoder->range_decoder, &decoder->bm_rep1[*state] ) == 0 )
            distance = decoder->rep1;
          else
            {
            if( Rd_decode_bit( decoder->range_decoder, &decoder->bm_rep2[*state] ) == 0 )
              distance = decoder->rep2;
            else { distance = decoder->rep3; decoder->rep3 = decoder->rep2; }
            decoder->rep2 = decoder->rep1;
            }
          decoder->rep1 = decoder->rep0;
          decoder->rep0 = distance;
          }
        else
          {
          if( Rd_decode_bit( decoder->range_decoder, &decoder->bm_len[*state][pos_state] ) == 0 )
            { St_set_short_rep( state ); len = 1; }
          }
        if( len == 0 )
          {
          St_set_rep( state );
          len = min_match_len + Led_decode( &decoder->rep_match_len_decoder, decoder->range_decoder, pos_state );
          }
        }
      else
        {
        int dis_slot;
        const unsigned int rep0_saved = decoder->rep0;
        len = min_match_len + Led_decode( &decoder->len_decoder, decoder->range_decoder, pos_state );
        dis_slot = Rd_decode_tree( decoder->range_decoder, decoder->bm_dis_slot[get_dis_state(len)], dis_slot_bits );
        if( dis_slot < start_dis_model ) decoder->rep0 = dis_slot;
        else
          {
          const int direct_bits = ( dis_slot >> 1 ) - 1;
          decoder->rep0 = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
          if( dis_slot < end_dis_model )
            decoder->rep0 += Rd_decode_tree_reversed( decoder->range_decoder, decoder->bm_dis + decoder->rep0 - dis_slot, direct_bits );
          else
            {
            decoder->rep0 += Rd_decode( decoder->range_decoder, direct_bits - dis_align_bits ) << dis_align_bits;
            decoder->rep0 += Rd_decode_tree_reversed( decoder->range_decoder, decoder->bm_align, dis_align_bits );
            if( decoder->rep0 == 0xFFFFFFFFU )		/* Marker found */
              {
              decoder->rep0 = rep0_saved;
              Rd_normalize( decoder->range_decoder );
              if( len == min_match_len )	/* End Of Stream marker */
                {
                if( Rd_available_bytes( decoder->range_decoder ) < Ft_versioned_size( decoder->member_version ) &&
                    !decoder->range_decoder->at_stream_end )
                  { decoder->verify_trailer_pending = true; return 0; }
                decoder->member_finished = true;
                if( LZd_verify_trailer( decoder ) ) return 0; else return 3;
                }
              if( len == min_match_len + 1 )	/* Sync Flush marker */
                {
                if( Rd_try_reload( decoder->range_decoder, true ) ) continue;
                else return 0;
                }
              return 4;
              }
            }
          }
        decoder->rep3 = decoder->rep2;
        decoder->rep2 = decoder->rep1; decoder->rep1 = rep0_saved;
        St_set_match( state );
        if( decoder->rep0 >= (unsigned int)decoder->dictionary_size ||
            ( decoder->rep0 >= (unsigned int)decoder->cb.put &&
              !decoder->partial_data_pos ) )
          return 1;
        }
      LZd_copy_block( decoder, decoder->rep0, len );
      }
    }
  return 2;
  }
