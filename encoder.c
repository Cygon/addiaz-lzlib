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

enum { max_num_trials = 1 << 12,
       price_shift = 6 };

static inline int price0( const Bit_model probability )
  { return get_price( probability ); }

static inline int price1( const Bit_model probability )
  { return get_price( bit_model_total-probability ); }

static inline int price_bit( const Bit_model bm, const int bit )
  { if( bit ) return price1( bm ); else return price0( bm ); }


static inline int price_symbol( const Bit_model bm[], int symbol,
                                const int num_bits )
  {
  int price = 0;
  symbol |= ( 1 << num_bits );
  while( symbol > 1 )
    {
    const int bit = symbol & 1;
    symbol >>= 1;
    price += price_bit( bm[symbol], bit );
    }
  return price;
  }


static inline int price_symbol_reversed( const Bit_model bm[], int symbol,
                                         const int num_bits )
  {
  int price = 0;
  int model = 1;
  int i;
  for( i = num_bits; i > 0; --i )
    {
    const int bit = symbol & 1;
    symbol >>= 1;
    price += price_bit( bm[model], bit );
    model = ( model << 1 ) | bit;
    }
  return price;
  }


static inline int price_matched( const Bit_model bm[], const int symbol,
                                 const int match_byte )
  {
  int price = 0;
  int model = 1;
  int i;

  for( i = 7; i >= 0; --i )
    {
    const int match_bit = ( match_byte >> i ) & 1;
    int bit = ( symbol >> i ) & 1;
    price += price_bit( bm[(match_bit<<8)+model+0x100], bit );
    model = ( model << 1 ) | bit;
    if( match_bit != bit )
      {
      while( --i >= 0 )
        {
        bit = ( symbol >> i ) & 1;
        price += price_bit( bm[model], bit );
        model = ( model << 1 ) | bit;
        }
      break;
      }
    }
  return price;
  }


enum { /* bytes to keep in buffer before dictionary */
       before_size = max_num_trials + 1,
       /* bytes to keep in buffer after pos */
       after_size = max_num_trials + max_match_len,
       num_prev_positions4 = 1 << 20,
       num_prev_positions3 = 1 << 18,
       num_prev_positions2 = 1 << 16,
       num_prev_positions = num_prev_positions4 + num_prev_positions3 +
                            num_prev_positions2 };

struct Matchfinder
  {
  long long partial_data_pos;
  uint8_t * buffer;		/* input buffer */
  int32_t * prev_positions;	/* last seen position of key */
  int32_t * prev_pos_tree;
  int dictionary_size;		/* bytes to keep in buffer before pos */
  int buffer_size;
  int pos;			/* current pos in buffer */
  int cyclic_pos;		/* current pos in dictionary */
  int stream_pos;		/* first byte not yet read from file */
  int pos_limit;		/* when reached, a new block must be read */
  int match_len_limit;
  int cycles;
  bool at_stream_end;		/* stream_pos shows real end of file */
  bool been_flushed;
  };


static int Mf_write_data( struct Matchfinder * const mf,
                          const uint8_t * const inbuf, const int size )
  {
  const int sz = min( mf->buffer_size - mf->stream_pos, size );
  if( !mf->at_stream_end && sz > 0 )
    {
    memcpy( mf->buffer + mf->stream_pos, inbuf, sz );
    mf->stream_pos += sz;
    }
  return sz;
  }


static bool Mf_normalize_pos( struct Matchfinder * const mf )
  {
  if( mf->pos > mf->stream_pos )
    { mf->pos = mf->stream_pos; return false; }
  if( !mf->at_stream_end )
    {
    int i;
    const int offset = mf->pos - mf->dictionary_size - before_size;
    const int size = mf->stream_pos - offset;
    memmove( mf->buffer, mf->buffer + offset, size );
    mf->partial_data_pos += offset;
    mf->pos -= offset;
    mf->stream_pos -= offset;
    for( i = 0; i < num_prev_positions; ++i )
      if( mf->prev_positions[i] >= 0 ) mf->prev_positions[i] -= offset;
    for( i = 0; i < 2 * mf->dictionary_size; ++i )
      if( mf->prev_pos_tree[i] >= 0 ) mf->prev_pos_tree[i] -= offset;
    }
  return true;
  }


static bool Mf_init( struct Matchfinder * const mf,
                     const int dict_size, const int len_limit )
  {
  int i;
  mf->partial_data_pos = 0;
  mf->dictionary_size = dict_size;
  mf->buffer_size = ( 2 * max( 65536, dict_size ) ) + before_size + after_size;
  mf->buffer = (uint8_t *)malloc( mf->buffer_size );
  mf->prev_positions = (int32_t *)malloc( num_prev_positions * sizeof (int32_t) );
  mf->prev_pos_tree = (int32_t *)malloc( 2 * dict_size * sizeof (int32_t) );
  mf->pos = 0;
  mf->cyclic_pos = 0;
  mf->stream_pos = 0;
  mf->pos_limit = mf->buffer_size - after_size;
  mf->match_len_limit = len_limit;
  mf->cycles = ( len_limit < max_match_len ) ? 16 + ( len_limit / 2 ) : 256;
  mf->at_stream_end = false;
  mf->been_flushed = false;

  if( !mf->buffer || !mf->prev_positions || !mf->prev_pos_tree )
    {
    if( mf->prev_pos_tree )
      { free( mf->prev_pos_tree ); mf->prev_pos_tree = 0; }
    if( mf->prev_positions )
      { free( mf->prev_positions ); mf->prev_positions = 0; }
    if( mf->buffer ) { free( mf->buffer ); mf->buffer = 0; }
    return false;
    }
  for( i = 0; i < num_prev_positions; ++i ) mf->prev_positions[i] = -1;
  return true;
  }


static inline void Mf_free( struct Matchfinder * const mf )
  {
  free( mf->prev_pos_tree ); mf->prev_pos_tree = 0;
  free( mf->prev_positions ); mf->prev_positions = 0;
  free( mf->buffer ); mf->buffer = 0;
  }

static inline uint8_t Mf_peek( const struct Matchfinder * const mf, const int i )
  { return mf->buffer[mf->pos+i]; }

static inline int Mf_available_bytes( const struct Matchfinder * const mf )
  { return mf->stream_pos - mf->pos; }

static inline long long Mf_data_position( const struct Matchfinder * const mf )
  { return mf->partial_data_pos + mf->pos; }

static inline bool Mf_finished( const struct Matchfinder * const mf )
  { return mf->at_stream_end && mf->pos >= mf->stream_pos; }

static inline const uint8_t * Mf_ptr_to_current_pos( const struct Matchfinder * const mf )
  { return mf->buffer + mf->pos; }

static inline void Mf_set_flushing( struct Matchfinder * const mf,
                                    const bool flushing )
  { mf->at_stream_end = flushing; }

static inline int Mf_free_bytes( const struct Matchfinder * const mf )
  { if( mf->at_stream_end ) return 0; return mf->buffer_size - mf->stream_pos; }

static inline bool Mf_enough_available_bytes( const struct Matchfinder * const mf )
  {
  return ( mf->pos + after_size <= mf->stream_pos ||
           ( mf->at_stream_end && mf->pos < mf->stream_pos ) );
  }

static inline bool Mf_dec_pos( struct Matchfinder * const mf,
                               const int ahead )
  {
  if( ahead < 0 || mf->pos < ahead ) return false;
  mf->pos -= ahead;
  mf->cyclic_pos -= ahead;
  if( mf->cyclic_pos < 0 ) mf->cyclic_pos += mf->dictionary_size;
  return true;
  }

static inline int Mf_true_match_len( const struct Matchfinder * const mf,
                                     const int index, const int distance,
                                     int len_limit )
  {
  const uint8_t * const data = mf->buffer + mf->pos + index;
  int i = 0;

  if( index + len_limit > Mf_available_bytes( mf ) )
    len_limit = Mf_available_bytes( mf ) - index;
  while( i < len_limit && data[i-distance] == data[i] ) ++i;
  return i;
  }

static inline bool Mf_move_pos( struct Matchfinder * const mf )
  {
  if( ++mf->cyclic_pos >= mf->dictionary_size ) mf->cyclic_pos = 0;
  if( ++mf->pos >= mf->pos_limit ) return Mf_normalize_pos( mf );
  return true;
  }


static void Mf_reset( struct Matchfinder * const mf )
  {
  int i;
  const int size = mf->stream_pos - mf->pos;
  if( size > 0 ) memmove( mf->buffer, mf->buffer + mf->pos, size );
  mf->partial_data_pos = 0;
  mf->stream_pos -= mf->pos;
  mf->pos = 0;
  mf->cyclic_pos = 0;
  mf->at_stream_end = false;
  mf->been_flushed = false;
  for( i = 0; i < num_prev_positions; ++i ) mf->prev_positions[i] = -1;
  }


static int Mf_longest_match_len( struct Matchfinder * const mf,
                                 int * const distances )
  {
  int32_t * ptr0 = mf->prev_pos_tree + ( mf->cyclic_pos << 1 );
  int32_t * ptr1 = ptr0 + 1;
  int32_t * newptr;
  const uint8_t * newdata;
  int len = 0, len0 = 0, len1 = 0;
  int maxlen = min_match_len - 1;
  const int min_pos = (mf->pos >= mf->dictionary_size) ?
                      (mf->pos - mf->dictionary_size + 1) : 0;
  const uint8_t * const data = mf->buffer + mf->pos;
  int count, delta, key2, key3, key4, newpos, tmp;
  int len_limit = mf->match_len_limit;

  if( len_limit > Mf_available_bytes( mf ) )
    {
    mf->been_flushed = true;
    len_limit = Mf_available_bytes( mf );
    if( len_limit < 4 ) { *ptr0 = *ptr1 = -1; return 0; }
    }

  key2 = num_prev_positions4 + num_prev_positions3 +
         ( ( (int)data[0] << 8 ) | data[1] );
  tmp = crc32[data[0]] ^ data[1] ^ ( (uint32_t)data[2] << 8 );
  key3 = num_prev_positions4 + (int)( tmp & ( num_prev_positions3 - 1 ) );
  key4 = (int)( ( tmp ^ ( crc32[data[3]] << 5 ) ) &
                ( num_prev_positions4 - 1 ) );

  if( distances )
    {
    int np = mf->prev_positions[key2];
    if( np >= min_pos )
      { distances[2] = mf->pos - np - 1; maxlen = 2; }
    else distances[2] = 0x7FFFFFFF;
    np = mf->prev_positions[key3];
    if( np >= min_pos && mf->buffer[np] == data[0] )
      { distances[3] = mf->pos - np - 1; maxlen = 3; }
    else distances[3] = 0x7FFFFFFF;
    distances[4] = 0x7FFFFFFF;
    }

  mf->prev_positions[key2] = mf->pos;
  mf->prev_positions[key3] = mf->pos;
  newpos = mf->prev_positions[key4];
  mf->prev_positions[key4] = mf->pos;

  for( count = mf->cycles; ; )
    {
    if( newpos < min_pos || --count < 0 ) { *ptr0 = *ptr1 = -1; break; }
    newdata = mf->buffer + newpos;
    if( mf->been_flushed ) len = 0;
    while( len < len_limit && newdata[len] == data[len] ) ++len;

    delta = mf->pos - newpos;
    if( distances ) while( maxlen < len ) distances[++maxlen] = delta - 1;

    newptr = mf->prev_pos_tree +
      ( ( mf->cyclic_pos - delta +
          ( ( mf->cyclic_pos >= delta ) ? 0 : mf->dictionary_size ) ) << 1 );

    if( len < len_limit )
      {
      if( newdata[len] < data[len] )
        {
        *ptr0 = newpos;
        ptr0 = newptr + 1;
        newpos = *ptr0;
        len0 = len; if( len1 < len ) len = len1;
        }
      else
        {
        *ptr1 = newpos;
        ptr1 = newptr;
        newpos = *ptr1;
        len1 = len; if( len0 < len ) len = len0;
        }
      }
    else
      {
      *ptr0 = newptr[0];
      *ptr1 = newptr[1];
      break;
      }
    }
  if( distances )
    {
    if( distances[3] > distances[4] ) distances[3] = distances[4];
    if( distances[2] > distances[3] ) distances[2] = distances[3];
    }
  return maxlen;
  }


enum { re_min_free_bytes = 2 * max_num_trials };

struct Range_encoder
  {
  struct Circular_buffer cb;
  uint64_t low;
  long long partial_member_pos;
  uint32_t range;
  int ff_count;
  uint8_t cache;
  };

static inline void Re_shift_low( struct Range_encoder * const renc )
  {
  const bool carry = ( renc->low > 0xFFFFFFFFU );
  if( carry || renc->low < 0xFF000000U )
    {
    Cb_put_byte( &renc->cb, renc->cache + carry );
    for( ; renc->ff_count > 0; --renc->ff_count )
      Cb_put_byte( &renc->cb, 0xFF + carry );
    renc->cache = renc->low >> 24;
    }
  else ++renc->ff_count;
  renc->low = ( renc->low & 0x00FFFFFFU ) << 8;
  }

static inline bool Re_init( struct Range_encoder * const renc )
  {
  if( !Cb_init( &renc->cb, 65536 + re_min_free_bytes ) ) return false;
  renc->low = 0;
  renc->partial_member_pos = 0;
  renc->range = 0xFFFFFFFFU;
  renc->ff_count = 0;
  renc->cache = 0;
  return true;
  }

static inline void Re_free( struct Range_encoder * const renc )
  { Cb_free( &renc->cb ); }

static inline long long Re_member_position( const struct Range_encoder * const renc )
  { return renc->partial_member_pos + Cb_used_bytes( &renc->cb ) + renc->ff_count; }

static inline bool Re_enough_free_bytes( const struct Range_encoder * const renc )
  { return Cb_free_bytes( &renc->cb ) >= re_min_free_bytes; }

static inline int Re_read_data( struct Range_encoder * const renc,
                                uint8_t * const out_buffer, const int out_size )
  {
  const int size = Cb_read_data( &renc->cb, out_buffer, out_size );
  if( size > 0 ) renc->partial_member_pos += size;
  return size;
  }

static inline void Re_flush( struct Range_encoder * const renc )
  {
  int i; for( i = 0; i < 5; ++i ) Re_shift_low( renc );
  renc->low = 0;
  renc->range = 0xFFFFFFFFU;
  renc->ff_count = 0;
  renc->cache = 0;
  }

static inline void Re_encode( struct Range_encoder * const renc,
                              const int symbol, const int num_bits )
  {
  int i;
  for( i = num_bits - 1; i >= 0; --i )
    {
    renc->range >>= 1;
    if( (symbol >> i) & 1 ) renc->low += renc->range;
    if( renc->range <= 0x00FFFFFFU )
      { renc->range <<= 8; Re_shift_low( renc ); }
    }
  }

static inline void Re_encode_bit( struct Range_encoder * const renc,
                                  Bit_model * const probability, const int bit )
  {
  const uint32_t bound = ( renc->range >> bit_model_total_bits ) * *probability;
  if( !bit )
    {
    renc->range = bound;
    *probability += (bit_model_total - *probability) >> bit_model_move_bits;
    }
  else
    {
    renc->low += bound;
    renc->range -= bound;
    *probability -= *probability >> bit_model_move_bits;
    }
  if( renc->range <= 0x00FFFFFFU )
    { renc->range <<= 8; Re_shift_low( renc ); }
  }

static inline void Re_encode_tree( struct Range_encoder * const renc,
                                   Bit_model bm[], const int symbol, const int num_bits )
  {
  int mask = ( 1 << ( num_bits - 1 ) );
  int model = 1;
  int i;
  for( i = num_bits; i > 0; --i, mask >>= 1 )
    {
    const int bit = ( symbol & mask );
    Re_encode_bit( renc, &bm[model], bit );
    model <<= 1;
    if( bit ) model |= 1;
    }
  }

static inline void Re_encode_tree_reversed( struct Range_encoder * const renc,
                                            Bit_model bm[], int symbol, const int num_bits )
  {
  int model = 1;
  int i;
  for( i = num_bits; i > 0; --i )
    {
    const int bit = symbol & 1;
    Re_encode_bit( renc, &bm[model], bit );
    model = ( model << 1 ) | bit;
    symbol >>= 1;
    }
  }

static inline void Re_encode_matched( struct Range_encoder * const renc,
                                      Bit_model bm[], int symbol, int match_byte )
  {
  int model = 1;
  int i;
  for( i = 7; i >= 0; --i )
    {
    const int match_bit = ( match_byte >> i ) & 1;
    int bit = ( symbol >> i ) & 1;
    Re_encode_bit( renc, &bm[(match_bit<<8)+model+0x100], bit );
    model = ( model << 1 ) | bit;
    if( match_bit != bit )
      {
      while( --i >= 0 )
        {
        bit = ( symbol >> i ) & 1;
        Re_encode_bit( renc, &bm[model], bit );
        model = ( model << 1 ) | bit;
        }
      break;
      }
    }
  }


struct Len_encoder
  {
  Bit_model choice1;
  Bit_model choice2;
  Bit_model bm_low[pos_states][len_low_symbols];
  Bit_model bm_mid[pos_states][len_mid_symbols];
  Bit_model bm_high[len_high_symbols];
  int prices[pos_states][max_len_symbols];
  int len_symbols;
  int counters[pos_states];
  };

static void Lee_update_prices( struct Len_encoder * const len_encoder,
                               const int pos_state )
  {
  int * const pps = len_encoder->prices[pos_state];
  int tmp = price0( len_encoder->choice1 );
  int len = 0;
  for( ; len < len_low_symbols && len < len_encoder->len_symbols; ++len )
    pps[len] = tmp +
               price_symbol( len_encoder->bm_low[pos_state], len, len_low_bits );
  tmp = price1( len_encoder->choice1 );
  for( ; len < len_low_symbols + len_mid_symbols && len < len_encoder->len_symbols; ++len )
    pps[len] = tmp + price0( len_encoder->choice2 ) +
               price_symbol( len_encoder->bm_mid[pos_state], len - len_low_symbols, len_mid_bits );
  for( ; len < len_encoder->len_symbols; ++len )
    /* using 4 slots per value makes "Lee_price" faster */
    len_encoder->prices[3][len] = len_encoder->prices[2][len] =
    len_encoder->prices[1][len] = len_encoder->prices[0][len] =
      tmp + price1( len_encoder->choice2 ) +
      price_symbol( len_encoder->bm_high, len - len_low_symbols - len_mid_symbols, len_high_bits );
  len_encoder->counters[pos_state] = len_encoder->len_symbols;
  }

static void Lee_init( struct Len_encoder * const len_encoder,
                      const int len_limit )
  {
  int i, j;
  Bm_init( &len_encoder->choice1 );
  Bm_init( &len_encoder->choice2 );
  for( i = 0; i < pos_states; ++i )
    for( j = 0; j < len_low_symbols; ++j )
      Bm_init( &len_encoder->bm_low[i][j] );
  for( i = 0; i < pos_states; ++i )
    for( j = 0; j < len_mid_symbols; ++j )
      Bm_init( &len_encoder->bm_mid[i][j] );
  for( i = 0; i < len_high_symbols; ++i )
    Bm_init( &len_encoder->bm_high[i] );
  len_encoder->len_symbols = len_limit + 1 - min_match_len;
  for( i = 0; i < pos_states; ++i ) Lee_update_prices( len_encoder, i );
  }


static void Lee_encode( struct Len_encoder * const len_encoder,
                        struct Range_encoder * const renc,
                        int symbol, const int pos_state )
  {
  symbol -= min_match_len;
  if( symbol < len_low_symbols )
    {
    Re_encode_bit( renc, &len_encoder->choice1, 0 );
    Re_encode_tree( renc, len_encoder->bm_low[pos_state], symbol, len_low_bits );
    }
  else
    {
    Re_encode_bit( renc, &len_encoder->choice1, 1 );
    if( symbol < len_low_symbols + len_mid_symbols )
      {
      Re_encode_bit( renc, &len_encoder->choice2, 0 );
      Re_encode_tree( renc, len_encoder->bm_mid[pos_state],
                      symbol - len_low_symbols, len_mid_bits );
      }
    else
      {
      Re_encode_bit( renc, &len_encoder->choice2, 1 );
      Re_encode_tree( renc, len_encoder->bm_high,
                      symbol - len_low_symbols - len_mid_symbols, len_high_bits );
      }
    }
  if( --len_encoder->counters[pos_state] <= 0 )
    Lee_update_prices( len_encoder, pos_state );
  }


static inline int Lee_price( const struct Len_encoder * const len_encoder,
                             const int symbol, const int pos_state )
  { return len_encoder->prices[pos_state][symbol - min_match_len]; }


struct Literal_encoder
  {
  Bit_model bm_literal[1<<literal_context_bits][0x300];
  };

static inline void Lie_init( struct Literal_encoder * const lienc )
  {
  int i, j;
  for( i = 0; i < 1<<literal_context_bits; ++i )
    for( j = 0; j < 0x300; ++j )
      Bm_init( &lienc->bm_literal[i][j] );
  }

static inline int Lie_state( const uint8_t prev_byte )
  { return ( prev_byte >> ( 8 - literal_context_bits ) ); }

static inline void Lie_encode( struct Literal_encoder * const lienc,
                               struct Range_encoder * const renc,
                               uint8_t prev_byte, uint8_t symbol )
  { Re_encode_tree( renc, lienc->bm_literal[Lie_state(prev_byte)], symbol, 8 ); }

static inline void Lie_encode_matched( struct Literal_encoder * const lienc,
                                       struct Range_encoder * const renc,
                                       uint8_t prev_byte, uint8_t symbol,
                                       uint8_t match_byte )
  { Re_encode_matched( renc, lienc->bm_literal[Lie_state(prev_byte)],
                       symbol, match_byte ); }

static inline int Lie_price_symbol( const struct Literal_encoder * const lienc,
                                    uint8_t prev_byte, uint8_t symbol )
  { return price_symbol( lienc->bm_literal[Lie_state(prev_byte)], symbol, 8 ); }

static inline int Lie_price_matched( const struct Literal_encoder * const lienc,
                                     uint8_t prev_byte, uint8_t symbol,
                                     uint8_t match_byte )
  { return price_matched( lienc->bm_literal[Lie_state(prev_byte)],
                          symbol, match_byte ); }


enum { infinite_price = 0x0FFFFFFF,
       max_marker_size = 16,
       num_rep_distances = 4 };		/* must be 4 */

struct Trial
  {
  State state;
  int dis;
  int prev_index;	/* index of prev trial in trials[] */
  int price;		/* dual use var; cumulative price, match length */
  int reps[num_rep_distances];
  };

static inline void Tr_update( struct Trial * const trial,
                              const int d, const int p_i, const int pr )
  {
  if( pr < trial->price )
    { trial->dis = d; trial->prev_index = p_i; trial->price = pr; }
  }


struct LZ_encoder
  {
  long long member_size_limit;
  int longest_match_found;
  uint32_t crc;

  Bit_model bm_match[states][pos_states];
  Bit_model bm_rep[states];
  Bit_model bm_rep0[states];
  Bit_model bm_rep1[states];
  Bit_model bm_rep2[states];
  Bit_model bm_len[states][pos_states];
  Bit_model bm_dis_slot[max_dis_states][1<<dis_slot_bits];
  Bit_model bm_dis[modeled_distances-end_dis_model+1];
  Bit_model bm_align[dis_align_size];

  struct Matchfinder * matchfinder;
  struct Range_encoder range_encoder;
  struct Len_encoder len_encoder;
  struct Len_encoder rep_match_len_encoder;
  struct Literal_encoder literal_encoder;

  int num_dis_slots;
  int rep_distances[num_rep_distances];
  int match_distances[max_match_len+1];
  struct Trial trials[max_num_trials];

  int dis_slot_prices[max_dis_states][2*max_dictionary_bits];
  int dis_prices[max_dis_states][modeled_distances];
  int align_prices[dis_align_size];
  int align_price_count;
  int fill_counter;
  State state;
  bool member_finished;
  };


static inline bool LZe_member_finished( const struct LZ_encoder * const encoder )
  {
  return ( encoder->member_finished &&
           !Cb_used_bytes( &encoder->range_encoder.cb ) );
  }


static void LZe_fill_align_prices( struct LZ_encoder * const encoder )
  {
  int i;
  for( i = 0; i < dis_align_size; ++i )
    encoder->align_prices[i] =
      price_symbol_reversed( encoder->bm_align, i, dis_align_bits );
  encoder->align_price_count = dis_align_size;
  }


static void LZe_fill_distance_prices( struct LZ_encoder * const encoder )
  {
  int dis, dis_state;
  for( dis = start_dis_model; dis < modeled_distances; ++dis )
    {
    const int dis_slot = dis_slots[dis];
    const int direct_bits = ( dis_slot >> 1 ) - 1;
    const int base = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
    const int price =
      price_symbol_reversed( encoder->bm_dis + base - dis_slot,
                             dis - base, direct_bits );
    for( dis_state = 0; dis_state < max_dis_states; ++dis_state )
      encoder->dis_prices[dis_state][dis] = price;
    }

  for( dis_state = 0; dis_state < max_dis_states; ++dis_state )
    {
    int * const dsp = encoder->dis_slot_prices[dis_state];
    int * const dp = encoder->dis_prices[dis_state];
    const Bit_model * const bmds = encoder->bm_dis_slot[dis_state];
    int slot = 0;
    for( ; slot < end_dis_model && slot < encoder->num_dis_slots; ++slot )
      dsp[slot] = price_symbol( bmds, slot, dis_slot_bits );
    for( ; slot < encoder->num_dis_slots; ++slot )
      dsp[slot] = price_symbol( bmds, slot, dis_slot_bits ) +
                  (((( slot >> 1 ) - 1 ) - dis_align_bits ) << price_shift );

    for( dis = 0; dis < start_dis_model; ++dis )
      dp[dis] = dsp[dis];
    for( ; dis < modeled_distances; ++dis )
      dp[dis] += dsp[dis_slots[dis]];
    }
  }


static bool LZe_init( struct LZ_encoder * const encoder,
                      struct Matchfinder * const mf,
                      const File_header header, const long long member_size )
  {
  int i, j;
  encoder->member_size_limit = member_size - Ft_size - max_marker_size;
  encoder->longest_match_found = 0;
  encoder->crc = 0xFFFFFFFFU;

  for( i = 0; i < states; ++i )
    {
    for( j = 0; j < pos_states; ++j )
      {
      Bm_init( &encoder->bm_match[i][j] );
      Bm_init( &encoder->bm_len[i][j] );
      }
    Bm_init( &encoder->bm_rep[i] );
    Bm_init( &encoder->bm_rep0[i] );
    Bm_init( &encoder->bm_rep1[i] );
    Bm_init( &encoder->bm_rep2[i] );
    }
  for( i = 0; i < max_dis_states; ++i )
    for( j = 0; j < 1<<dis_slot_bits; ++j )
      Bm_init( &encoder->bm_dis_slot[i][j] );
  for( i = 0; i < modeled_distances-end_dis_model+1; ++i )
    Bm_init( &encoder->bm_dis[i] );
  for( i = 0; i < dis_align_size; ++i )
    Bm_init( &encoder->bm_align[i] );

  encoder->matchfinder = mf;
  if( !Re_init( &encoder->range_encoder ) ) return false;
  Lee_init( &encoder->len_encoder, encoder->matchfinder->match_len_limit );
  Lee_init( &encoder->rep_match_len_encoder, encoder->matchfinder->match_len_limit );
  Lie_init( &encoder->literal_encoder );
  encoder->num_dis_slots =
    2 * real_bits( encoder->matchfinder->dictionary_size - 1 );
  for( i = 0; i < num_rep_distances; ++i ) encoder->rep_distances[i] = 0;
  encoder->fill_counter = 0;
  encoder->state = 0;
  encoder->member_finished = false;

  LZe_fill_align_prices( encoder );

  for( i = 0; i < Fh_size; ++i )
    Cb_put_byte( &encoder->range_encoder.cb, header[i] );
  return true;
  }


static inline void LZe_free( struct LZ_encoder * const encoder )
  {
  Re_free( &encoder->range_encoder );
  }

static inline uint32_t LZe_crc( const struct LZ_encoder * const encoder )
  { return encoder->crc ^ 0xFFFFFFFFU; }

       /* move-to-front dis in/into reps */
static inline void LZe_mtf_reps( const int dis, int reps[num_rep_distances] )
  {
  int i;
  if( dis >= num_rep_distances )
    {
    for( i = num_rep_distances - 1; i > 0; --i ) reps[i] = reps[i-1];
    reps[0] = dis - num_rep_distances;
    }
  else if( dis > 0 )
    {
    const int distance = reps[dis];
    for( i = dis; i > 0; --i ) reps[i] = reps[i-1];
    reps[0] = distance;
    }
  }

static inline int LZe_price_rep_len1( const struct LZ_encoder * const encoder,
                                      const State state, const int pos_state )
  {
  return price0( encoder->bm_rep0[state] ) +
         price0( encoder->bm_len[state][pos_state] );
  }

static inline int LZe_price_rep( const struct LZ_encoder * const encoder,
                                 const int rep,
                                 const State state, const int pos_state )
  {
  int price;
  if( rep == 0 ) return price0( encoder->bm_rep0[state] ) +
                        price1( encoder->bm_len[state][pos_state] );
  price = price1( encoder->bm_rep0[state] );
  if( rep == 1 )
    price += price0( encoder->bm_rep1[state] );
  else
    {
    price += price1( encoder->bm_rep1[state] );
    price += price_bit( encoder->bm_rep2[state], rep - 2 );
    }
  return price;
  }

static inline int LZe_price_dis( const struct LZ_encoder * const encoder,
                                 const int dis, const int dis_state )
  {
  if( dis < modeled_distances )
    return encoder->dis_prices[dis_state][dis];
  else
    return encoder->dis_slot_prices[dis_state][get_slot( dis )] +
           encoder->align_prices[dis & (dis_align_size - 1)];
  }

static inline int LZe_price_pair( const struct LZ_encoder * const encoder,
                                  const int dis, const int len,
                                  const int pos_state )
  {
  if( len <= min_match_len && dis >= modeled_distances )
    return infinite_price;
  return Lee_price( &encoder->len_encoder, len, pos_state ) +
         LZe_price_dis( encoder, dis, get_dis_state( len ) );
  }

static inline void LZe_encode_pair( struct LZ_encoder * const encoder,
                                    const uint32_t dis, const int len,
                                    const int pos_state )
  {
  const int dis_slot = get_slot( dis );
  Lee_encode( &encoder->len_encoder, &encoder->range_encoder, len, pos_state );
  Re_encode_tree( &encoder->range_encoder,
                  encoder->bm_dis_slot[get_dis_state(len)],
                  dis_slot, dis_slot_bits );

  if( dis_slot >= start_dis_model )
    {
    const int direct_bits = ( dis_slot >> 1 ) - 1;
    const uint32_t base = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
    const uint32_t direct_dis = dis - base;

    if( dis_slot < end_dis_model )
      Re_encode_tree_reversed( &encoder->range_encoder,
                               encoder->bm_dis + base - dis_slot,
                               direct_dis, direct_bits );
    else
      {
      Re_encode( &encoder->range_encoder, direct_dis >> dis_align_bits,
                 direct_bits - dis_align_bits );
      Re_encode_tree_reversed( &encoder->range_encoder, encoder->bm_align,
                               direct_dis, dis_align_bits );
      if( --encoder->align_price_count <= 0 ) LZe_fill_align_prices( encoder );
      }
    }
  }

     /* End Of Stream mark => (dis == 0xFFFFFFFFU, len == min_match_len) */
static bool LZe_full_flush( struct LZ_encoder * const encoder, const State state )
  {
  int i;
  const int pos_state = Mf_data_position( encoder->matchfinder ) & pos_state_mask;
  File_trailer trailer;
  if( encoder->member_finished ||
      Cb_free_bytes( &encoder->range_encoder.cb ) < Ft_size + max_marker_size )
    return false;
  Re_encode_bit( &encoder->range_encoder, &encoder->bm_match[state][pos_state], 1 );
  Re_encode_bit( &encoder->range_encoder, &encoder->bm_rep[state], 0 );
  LZe_encode_pair( encoder, 0xFFFFFFFFU, min_match_len, pos_state );
  Re_flush( &encoder->range_encoder );
  Ft_set_data_crc( trailer, LZe_crc( encoder ) );
  Ft_set_data_size( trailer, Mf_data_position( encoder->matchfinder ) );
  Ft_set_member_size( trailer, Re_member_position( &encoder->range_encoder ) +
                               Ft_size );
  for( i = 0; i < Ft_size; ++i )
    Cb_put_byte( &encoder->range_encoder.cb, trailer[i] );
  return true;
  }


     /* Sync Flush mark => (dis == 0xFFFFFFFFU, len == min_match_len + 1) */
static bool LZe_sync_flush( struct LZ_encoder * const encoder )
  {
  const int pos_state = Mf_data_position( encoder->matchfinder ) & pos_state_mask;
  const State state = encoder->state;
  if( encoder->member_finished ||
      Cb_free_bytes( &encoder->range_encoder.cb ) < max_marker_size )
    return false;
  Re_encode_bit( &encoder->range_encoder, &encoder->bm_match[state][pos_state], 1 );
  Re_encode_bit( &encoder->range_encoder, &encoder->bm_rep[state], 0 );
  LZe_encode_pair( encoder, 0xFFFFFFFFU, min_match_len + 1, pos_state );
  Re_flush( &encoder->range_encoder );
  return true;
  }


static inline int LZe_read_match_distances( struct LZ_encoder * const encoder )
  {
  int len = Mf_longest_match_len( encoder->matchfinder,
                                  encoder->match_distances );
  if( len == encoder->matchfinder->match_len_limit && len < max_match_len )
    len += Mf_true_match_len( encoder->matchfinder, len,
                              encoder->match_distances[len] + 1,
                              max_match_len - len );
  return len;
  }

static inline bool LZe_move_pos( struct LZ_encoder * const encoder, int n )
  {
  if( --n >= 0 && !Mf_move_pos( encoder->matchfinder ) ) return false;
  while( --n >= 0 )
    {
    Mf_longest_match_len( encoder->matchfinder, 0 );
    if( !Mf_move_pos( encoder->matchfinder ) ) return false;
    }
  return true;
  }

static inline void LZe_backward( struct LZ_encoder * const encoder, int cur )
  {
  int * const dis = &encoder->trials[cur].dis;
  while( cur > 0 )
    {
    const int prev_index = encoder->trials[cur].prev_index;
    struct Trial * const prev_trial = &encoder->trials[prev_index];
    prev_trial->price = cur - prev_index;			/* len */
    cur = *dis; *dis = prev_trial->dis; prev_trial->dis = cur;
    cur = prev_index;
    }
  }


/* Return value == number of bytes advanced (ahead).
   trials[0]..trials[retval-1] contain the steps to encode.
   ( trials[0].dis == -1 && trials[0].price == 1 ) means literal.
*/
static int LZe_sequence_optimizer( struct LZ_encoder * const encoder,
                                   const int reps[num_rep_distances],
                                   const State state )
  {
  int main_len, i, rep, cur = 0, num_trials;
  int replens[num_rep_distances];
  int rep_index = 0;

  if( encoder->longest_match_found > 0 )	/* from previous call */
    {
    main_len = encoder->longest_match_found;
    encoder->longest_match_found = 0;
    }
  else main_len = LZe_read_match_distances( encoder );

  for( i = 0; i < num_rep_distances; ++i )
    {
    replens[i] =
      Mf_true_match_len( encoder->matchfinder, 0, reps[i] + 1, max_match_len );
    if( replens[i] > replens[rep_index] ) rep_index = i;
    }
  if( replens[rep_index] >= encoder->matchfinder->match_len_limit )
    {
    encoder->trials[0].dis = rep_index;
    encoder->trials[0].price = replens[rep_index];
    if( !LZe_move_pos( encoder, replens[rep_index] ) ) return 0;
    return replens[rep_index];
    }

  if( main_len >= encoder->matchfinder->match_len_limit )
    {
    encoder->trials[0].dis =
      encoder->match_distances[encoder->matchfinder->match_len_limit] +
      num_rep_distances;
    encoder->trials[0].price = main_len;
    if( !LZe_move_pos( encoder, main_len ) ) return 0;
    return main_len;
    }

  {
  const int pos_state = Mf_data_position( encoder->matchfinder ) & pos_state_mask;
  const int match_price = price1( encoder->bm_match[state][pos_state] );
  const int rep_match_price = match_price + price1( encoder->bm_rep[state] );
  const uint8_t prev_byte = Mf_peek( encoder->matchfinder, -1 );
  const uint8_t cur_byte = Mf_peek( encoder->matchfinder, 0 );
  const uint8_t match_byte = Mf_peek( encoder->matchfinder, -reps[0]-1 );

  encoder->trials[0].state = state;
  for( i = 0; i < num_rep_distances; ++i )
    encoder->trials[0].reps[i] = reps[i];
  encoder->trials[1].dis = -1;
  encoder->trials[1].prev_index = 0;
  encoder->trials[1].price = price0( encoder->bm_match[state][pos_state] );
  if( St_is_char( state ) )
    encoder->trials[1].price +=
      Lie_price_symbol( &encoder->literal_encoder, prev_byte, cur_byte );
  else
    encoder->trials[1].price +=
      Lie_price_matched( &encoder->literal_encoder, prev_byte, cur_byte, match_byte );

  if( match_byte == cur_byte )
    Tr_update( &encoder->trials[1], 0, 0, rep_match_price +
               LZe_price_rep_len1( encoder, state, pos_state ) );

  if( main_len < min_match_len )
    {
    encoder->trials[0].dis = encoder->trials[1].dis;
    encoder->trials[0].price = 1;
    if( !Mf_move_pos( encoder->matchfinder ) ) return 0;
    return 1;
    }

  if( main_len <= replens[rep_index] )
    {
    int len;
    main_len = replens[rep_index];
    for( len = min_match_len; len <= main_len; ++len )
      encoder->trials[len].price = infinite_price;
    }
  else
    {
    int len;
    const int normal_match_price = match_price + price0( encoder->bm_rep[state] );
    for( len = min_match_len; len <= main_len; ++len )
      {
      encoder->trials[len].dis = encoder->match_distances[len] + num_rep_distances;
      encoder->trials[len].prev_index = 0;
      encoder->trials[len].price = normal_match_price +
        LZe_price_pair( encoder, encoder->match_distances[len], len, pos_state );
      }
    }

  for( rep = 0; rep < num_rep_distances; ++rep )
    {
    const int price = rep_match_price +
                      LZe_price_rep( encoder, rep, state, pos_state );
    int len;
    for( len = min_match_len; len <= replens[rep]; ++len )
      Tr_update( &encoder->trials[len], rep, 0, price +
                 Lee_price( &encoder->rep_match_len_encoder, len, pos_state ) );
    }
  }

  num_trials = main_len;
  if( !Mf_move_pos( encoder->matchfinder ) ) return 0;

  while( true )
    {
    struct Trial *cur_trial, *next_trial;
    int newlen, pos_state, prev_index, len_limit;
    int next_price, match_price, rep_match_price;
    uint8_t prev_byte, cur_byte, match_byte;

    if( ++cur >= num_trials )		/* no more initialized trials */
      {
      LZe_backward( encoder, cur );
      return cur;
      }
    newlen = LZe_read_match_distances( encoder );
    if( newlen >= encoder->matchfinder->match_len_limit )
      {
      encoder->longest_match_found = newlen;
      LZe_backward( encoder, cur );
      return cur;
      }

    cur_trial = &encoder->trials[cur];
    prev_index = cur_trial->prev_index;

    cur_trial->state = encoder->trials[prev_index].state;

    for( i = 0; i < num_rep_distances; ++i )
      cur_trial->reps[i] = encoder->trials[prev_index].reps[i];

    if( prev_index == cur - 1 )
      {
      if( cur_trial->dis == 0 ) St_set_short_rep( &cur_trial->state );
      else St_set_char( &cur_trial->state );
      }
    else
      {
      if( cur_trial->dis < num_rep_distances ) St_set_rep( &cur_trial->state );
      else St_set_match( &cur_trial->state );
      LZe_mtf_reps( cur_trial->dis, cur_trial->reps );
      }

    pos_state = Mf_data_position( encoder->matchfinder ) & pos_state_mask;
    prev_byte = Mf_peek( encoder->matchfinder, -1 );
    cur_byte = Mf_peek( encoder->matchfinder, 0 );
    match_byte = Mf_peek( encoder->matchfinder, -cur_trial->reps[0]-1 );

    next_price = cur_trial->price +
                 price0( encoder->bm_match[cur_trial->state][pos_state] );
    if( St_is_char( cur_trial->state ) )
      next_price += Lie_price_symbol( &encoder->literal_encoder,
                                      prev_byte, cur_byte );
    else
      next_price += Lie_price_matched( &encoder->literal_encoder,
                                       prev_byte, cur_byte, match_byte );
    if( !Mf_move_pos( encoder->matchfinder ) ) return 0;

    next_trial = &encoder->trials[cur+1];

    Tr_update( next_trial, -1, cur, next_price );

    match_price = cur_trial->price + price1( encoder->bm_match[cur_trial->state][pos_state] );
    rep_match_price = match_price + price1( encoder->bm_rep[cur_trial->state] );

    if( match_byte == cur_byte && next_trial->dis != 0 )
      Tr_update( next_trial, 0, cur, rep_match_price +
                 LZe_price_rep_len1( encoder, cur_trial->state, pos_state ) );

    len_limit = min( min( max_num_trials - 1 - cur,
                          Mf_available_bytes( encoder->matchfinder ) ),
                     encoder->matchfinder->match_len_limit );
    if( len_limit < min_match_len ) continue;

    for( rep = 0; rep < num_rep_distances; ++rep )
      {
      const int dis = cur_trial->reps[rep] + 1;
      int len = 0;
      const uint8_t * const data = Mf_ptr_to_current_pos( encoder->matchfinder ) - 1;
      while( len < len_limit && data[len] == data[len-dis] ) ++len;
      if( len >= min_match_len )
        {
        const int price = rep_match_price +
                          LZe_price_rep( encoder, rep, cur_trial->state, pos_state );
        while( num_trials < cur + len )
          encoder->trials[++num_trials].price = infinite_price;
        for( ; len >= min_match_len; --len )
          Tr_update( &encoder->trials[cur+len], rep, cur, price +
                     Lee_price( &encoder->rep_match_len_encoder, len, pos_state ) );
        }
      }

    if( newlen <= len_limit &&
        ( newlen > min_match_len ||
          ( newlen == min_match_len &&
            encoder->match_distances[min_match_len] < modeled_distances ) ) )
      {
      const int normal_match_price = match_price +
                                     price0( encoder->bm_rep[cur_trial->state] );
      int len;
      int dis = encoder->match_distances[min_match_len];
      int dis_state = get_dis_state( min_match_len );
      int dis_price = infinite_price;

      while( num_trials < cur + newlen )
        encoder->trials[++num_trials].price = infinite_price;

      if( dis < modeled_distances )
        Tr_update( &encoder->trials[cur+min_match_len], dis + num_rep_distances,
                   cur, normal_match_price + encoder->dis_prices[dis_state][dis] +
                   Lee_price( &encoder->len_encoder, min_match_len, pos_state ) );

      for( len = min_match_len + 1; len <= newlen; ++len )
        {
        if( dis != encoder->match_distances[len] ||
            dis_state < max_dis_states - 1 )
          {
          dis = encoder->match_distances[len];
          dis_state = get_dis_state( len );
          dis_price = LZe_price_dis( encoder, dis, dis_state );
          }
        Tr_update( &encoder->trials[cur+len], dis + num_rep_distances, cur,
                   normal_match_price + dis_price +
                   Lee_price( &encoder->len_encoder, len, pos_state ) );
        }
      }
    }
  }


static bool LZe_encode_member( struct LZ_encoder * const encoder,
                               const bool finish )
  {
  const int fill_count =
    ( encoder->matchfinder->match_len_limit > 12 ) ? 512 : 2048;
  int ahead, i;
  State * const state = &encoder->state;
  if( encoder->member_finished ) return true;
  if( Re_member_position( &encoder->range_encoder ) >= encoder->member_size_limit )
    {
    if( LZe_full_flush( encoder, *state ) ) encoder->member_finished = true;
    return true;
    }

  /* encode first byte */
  if( Mf_data_position( encoder->matchfinder ) == 0 &&
      !Mf_finished( encoder->matchfinder ) )
    {
    if( Mf_available_bytes( encoder->matchfinder ) < max_match_len &&
        !encoder->matchfinder->at_stream_end )
      return true;
    const uint8_t prev_byte = 0;
    const uint8_t cur_byte = Mf_peek( encoder->matchfinder, 0 );
    Re_encode_bit( &encoder->range_encoder, &encoder->bm_match[*state][0], 0 );
    Lie_encode( &encoder->literal_encoder, &encoder->range_encoder, prev_byte, cur_byte );
    CRC32_update_byte( &encoder->crc, cur_byte );
    Mf_longest_match_len( encoder->matchfinder, 0 );
    if( !Mf_move_pos( encoder->matchfinder ) ) return false;
    }

  while( !Mf_finished( encoder->matchfinder ) )
    {
    if( !Mf_enough_available_bytes( encoder->matchfinder ) ||
        !Re_enough_free_bytes( &encoder->range_encoder ) ) return true;
    if( encoder->fill_counter <= 0 )
      { LZe_fill_distance_prices( encoder ); encoder->fill_counter = fill_count; }

    ahead = LZe_sequence_optimizer( encoder, encoder->rep_distances, *state );
    if( ahead <= 0 ) return false;		/* can't happen */
    encoder->fill_counter -= ahead;

    for( i = 0; ; )
      {
      const int pos_state =
        ( Mf_data_position( encoder->matchfinder ) - ahead ) & pos_state_mask;
      const int dis = encoder->trials[i].dis;
      const int len = encoder->trials[i].price;

      bool bit = ( dis < 0 && len == 1 );
      Re_encode_bit( &encoder->range_encoder,
                     &encoder->bm_match[*state][pos_state], !bit );
      if( bit )					/* literal byte */
        {
        const uint8_t prev_byte = Mf_peek( encoder->matchfinder, -ahead-1 );
        const uint8_t cur_byte = Mf_peek( encoder->matchfinder, -ahead );
        CRC32_update_byte( &encoder->crc, cur_byte );
        if( St_is_char( *state ) )
          Lie_encode( &encoder->literal_encoder, &encoder->range_encoder,
                      prev_byte, cur_byte );
        else
          {
          const uint8_t match_byte =
            Mf_peek( encoder->matchfinder, -ahead-encoder->rep_distances[0]-1 );
          Lie_encode_matched( &encoder->literal_encoder, &encoder->range_encoder,
                              prev_byte, cur_byte, match_byte );
          }
        St_set_char( state );
        }
      else				/* match or repeated match */
        {
        CRC32_update_buf( &encoder->crc, Mf_ptr_to_current_pos( encoder->matchfinder ) - ahead, len );
        LZe_mtf_reps( dis, encoder->rep_distances );
        bit = ( dis < num_rep_distances );
        Re_encode_bit( &encoder->range_encoder, &encoder->bm_rep[*state], bit );
        if( bit )
          {
          bit = ( dis == 0 );
          Re_encode_bit( &encoder->range_encoder, &encoder->bm_rep0[*state], !bit );
          if( bit )
            Re_encode_bit( &encoder->range_encoder, &encoder->bm_len[*state][pos_state], len > 1 );
          else
            {
            Re_encode_bit( &encoder->range_encoder, &encoder->bm_rep1[*state], dis > 1 );
            if( dis > 1 )
              Re_encode_bit( &encoder->range_encoder, &encoder->bm_rep2[*state], dis > 2 );
            }
          if( len == 1 ) St_set_short_rep( state );
          else
            {
            Lee_encode( &encoder->rep_match_len_encoder, &encoder->range_encoder, len, pos_state );
            St_set_rep( state );
            }
          }
        else
          {
          LZe_encode_pair( encoder, dis - num_rep_distances, len, pos_state );
          St_set_match( state );
          }
        }
      ahead -= len; i += len;
      if( Re_member_position( &encoder->range_encoder ) >= encoder->member_size_limit )
        {
        if( !Mf_dec_pos( encoder->matchfinder, ahead ) ) return false;
        if( LZe_full_flush( encoder, *state ) ) encoder->member_finished = true;
        return true;
        }
      if( ahead <= 0 ) break;
      }
    }
  if( finish && LZe_full_flush( encoder, *state ) )
    encoder->member_finished = true;
  return true;
  }
