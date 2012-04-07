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

#ifndef max
  #define max(x,y) ((x) >= (y) ? (x) : (y))
#endif
#ifndef min
  #define min(x,y) ((x) <= (y) ? (x) : (y))
#endif

typedef uint8_t State;

enum { states = 12 };

static inline bool St_is_char( const State st ) { return st < 7; }

static inline void St_set_char( State * const st )
  {
  static const uint8_t next[states] =
    { 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 4, 5 };
  *st = next[*st];
  }

static inline void St_set_match( State * const st )
  {
  static const uint8_t next[states] =
    { 7, 7, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10 };
  *st = next[*st];
  }

static inline void St_set_rep( State * const st )
  {
  static const uint8_t next[states] =
    { 8, 8, 8, 8, 8, 8, 8, 11, 11, 11, 11, 11 };
  *st = next[*st];
  }

static inline void St_set_short_rep( State * const st )
  {
  static const uint8_t next[states] =
    { 9, 9, 9, 9, 9, 9, 9, 11, 11, 11, 11, 11 };
  *st = next[*st];
  }


enum {
  min_dictionary_bits = 12,
  min_dictionary_size = 1 << min_dictionary_bits,
  max_dictionary_bits = 29,
  max_dictionary_size = 1 << max_dictionary_bits,
  literal_context_bits = 3,
  pos_state_bits = 2,
  pos_states = 1 << pos_state_bits,
  pos_state_mask = pos_states - 1,

  dis_slot_bits = 6,
  start_dis_model = 4,
  end_dis_model = 14,
  modeled_distances = 1 << (end_dis_model / 2),
  dis_align_bits = 4,
  dis_align_size = 1 << dis_align_bits,

  len_low_bits = 3,
  len_mid_bits = 3,
  len_high_bits = 8,
  len_low_symbols = 1 << len_low_bits,
  len_mid_symbols = 1 << len_mid_bits,
  len_high_symbols = 1 << len_high_bits,
  max_len_symbols = len_low_symbols + len_mid_symbols + len_high_symbols,

  min_match_len = 2,					/* must be 2 */
  max_match_len = min_match_len + max_len_symbols - 1,	/* 273 */
  min_match_len_limit = 5,

  max_dis_states = 4 };

static inline int get_dis_state( int len )
  {
  len -= min_match_len;
  if( len >= max_dis_states ) len = max_dis_states - 1;
  return len;
  }


enum { bit_model_move_bits = 5,
       bit_model_total_bits = 11,
       bit_model_total = 1 << bit_model_total_bits };

typedef unsigned int Bit_model;

static inline void Bm_init( Bit_model * const probability )
  { *probability = bit_model_total / 2; }


static inline int real_bits( const unsigned int value )
  {
  int bits = 0, i = 1;
  unsigned int mask = 1;
  for( ; mask > 0; ++i, mask <<= 1 ) if( value & mask ) bits = i;
  return bits;
  }


static const uint8_t magic_string[4] = { 'L', 'Z', 'I', 'P' };

typedef uint8_t File_header[6];		/* 0-3 magic bytes */
					/*   4 version */
					/*   5 coded_dict_size */
enum { Fh_size = 6 };

static inline void Fh_set_magic( File_header data )
  { memcpy( data, magic_string, 4 ); data[4] = 1; }

static inline bool Fh_verify_magic( const File_header data )
  { return ( memcmp( data, magic_string, 4 ) == 0 ); }

static inline uint8_t Fh_version( const File_header data )
  { return data[4]; }

static inline bool Fh_verify_version( const File_header data )
  { return ( data[4] <= 1 ); }

static inline int Fh_get_dictionary_size( const File_header data )
  {
  int sz = ( 1 << ( data[5] & 0x1F ) );
  if( sz > min_dictionary_size && sz <= max_dictionary_size )
    sz -= ( sz / 16 ) * ( ( data[5] >> 5 ) & 0x07 );
  return sz;
  }

static inline bool Fh_set_dictionary_size( File_header data, const int sz )
  {
  if( sz >= min_dictionary_size && sz <= max_dictionary_size )
    {
    data[5] = real_bits( sz - 1 );
    if( sz > min_dictionary_size )
      {
      const int base_size = 1 << data[5];
      const int wedge = base_size / 16;
      int i;
      for( i = 7; i >= 1; --i )
        if( base_size - ( i * wedge ) >= sz )
          { data[5] |= ( i << 5 ); break; }
      }
    return true;
    }
  return false;
  }

static inline bool Fh_verify( const File_header data )
  {
  return ( Fh_verify_magic( data ) && Fh_verify_version( data ) &&
           Fh_get_dictionary_size( data ) >= min_dictionary_size &&
           Fh_get_dictionary_size( data ) <= max_dictionary_size );
  }


typedef uint8_t File_trailer[20];
			/*  0-3  CRC32 of the uncompressed data */
			/*  4-11 size of the uncompressed data */
			/* 12-19 member size including header and trailer */

enum { Ft_size = 20 };

static inline int Ft_versioned_size( const int version )
  { return ( ( version >= 1 ) ? 20 : 12 ); }

static inline uint32_t Ft_get_data_crc( const File_trailer data )
  {
  uint32_t tmp = 0;
  int i;
  for( i = 3; i >= 0; --i ) { tmp <<= 8; tmp += data[i]; }
  return tmp;
  }

static inline void Ft_set_data_crc( File_trailer data, uint32_t crc )
  {
  int i;
  for( i = 0; i <= 3; ++i ) { data[i] = (uint8_t)crc; crc >>= 8; }
  }

static inline long long Ft_get_data_size( const File_trailer data )
  {
  long long tmp = 0;
  int i;
  for( i = 11; i >= 4; --i ) { tmp <<= 8; tmp += data[i]; }
  return tmp;
  }

static inline void Ft_set_data_size( File_trailer data, long long sz )
  {
  int i;
  for( i = 4; i <= 11; ++i ) { data[i] = (uint8_t)sz; sz >>= 8; }
  }

static inline long long Ft_get_member_size( const File_trailer data )
  {
  long long tmp = 0;
  int i;
  for( i = 19; i >= 12; --i ) { tmp <<= 8; tmp += data[i]; }
  return tmp;
  }

static inline void Ft_set_member_size( File_trailer data, long long sz )
  {
  int i;
  for( i = 12; i <= 19; ++i ) { data[i] = (uint8_t)sz; sz >>= 8; }
  }
