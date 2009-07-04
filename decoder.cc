/*  Lzlib - A compression library for lzip files
    Copyright (C) 2009 Antonio Diaz Diaz.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

#define _FILE_OFFSET_BITS 64

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <stdint.h>

#include "lzlib.h"
#include "lzip.h"
#include "decoder.h"


const CRC32 crc32;


int Circular_buffer::read_data( uint8_t * const out_buffer, const int out_size ) throw()
  {
  int size = 0;
  if( get > put )
    {
    size = std::min( buffer_size - get, out_size );
    if( size > 0 )
      {
      std::memcpy( out_buffer, buffer + get, size );
      get += size;
      if( get >= buffer_size ) get = 0;
      }
    }
  if( get < put )
    {
    const int size2 = std::min( put - get, out_size - size );
    if( size2 > 0 )
      {
      std::memcpy( out_buffer + size, buffer + get, size2 );
      get += size2;
      size += size2;
      }
    }
  return size;
  }


int Circular_buffer::write_data( uint8_t * const in_buffer, const int in_size ) throw()
  {
  int size = 0;
  if( put >= get )
    {
    size = std::min( buffer_size - put - (get == 0), in_size );
    if( size > 0 )
      {
      std::memcpy( buffer + put, in_buffer, size );
      put += size;
      if( put >= buffer_size ) put = 0;
      }
    }
  if( put < get )
    {
    const int size2 = std::min( get - put - 1, in_size - size );
    if( size2 > 0 )
      {
      std::memcpy( buffer + put, in_buffer + size, size2 );
      put += size2;
      size += size2;
      }
    }
  return size;
  }


bool LZ_decoder::verify_trailer()
  {
  bool error = false;
  File_trailer trailer;
  const int trailer_size = trailer.size( format_version );
  for( int i = 0; i < trailer_size && !error; ++i )
    {
    if( !range_decoder.finished() )
      ((uint8_t *)&trailer)[i] = range_decoder.get_byte();
    else error = true;
    }
  if( format_version == 0 ) trailer.member_size( member_position() );
  if( trailer.data_crc() != crc() ) error = true;
  if( trailer.data_size() != data_position() ) error = true;
  if( trailer.member_size() != member_position() ) error = true;
  return !error;
  }


    // Return value: 0 = OK, 1 = decoder error, 2 = unexpected EOF,
    //               3 = trailer error, 4 = unknown marker found.
int LZ_decoder::decode_member()
  {
  if( member_finished_ ) return 0;
  if( !range_decoder.try_reload() ) return 0;
  while( true )
    {
    if( range_decoder.finished() ) return 2;
    if( !range_decoder.enough_available_bytes() || !enough_free_bytes() )
      return 0;
    const int pos_state = data_position() & pos_state_mask;
    if( range_decoder.decode_bit( bm_match[state()][pos_state] ) == 0 )
      {
      const uint8_t prev_byte = get_byte( 0 );
      if( state.is_char() )
        put_byte( literal_decoder.decode( range_decoder, prev_byte ) );
      else
        put_byte( literal_decoder.decode_matched( range_decoder, prev_byte,
                                                  get_byte( rep0 ) ) );
      state.set_char();
      }
    else
      {
      int len;
      if( range_decoder.decode_bit( bm_rep[state()] ) == 1 )
        {
        len = 0;
        if( range_decoder.decode_bit( bm_rep0[state()] ) == 0 )
          {
          if( range_decoder.decode_bit( bm_len[state()][pos_state] ) == 0 )
            { len = 1; state.set_short_rep(); }
          }
        else
          {
          unsigned int distance;
          if( range_decoder.decode_bit( bm_rep1[state()] ) == 0 )
            distance = rep1;
          else
            {
            if( range_decoder.decode_bit( bm_rep2[state()] ) == 0 )
              distance = rep2;
            else { distance = rep3; rep3 = rep2; }
            rep2 = rep1;
            }
          rep1 = rep0;
          rep0 = distance;
          }
        if( len == 0 )
          {
          len = min_match_len + rep_match_len_decoder.decode( range_decoder, pos_state );
          state.set_rep();
          }
        }
      else
        {
        unsigned int rep0_saved = rep0;
        len = min_match_len + len_decoder.decode( range_decoder, pos_state );
        const int dis_slot = range_decoder.decode_tree( bm_dis_slot[get_dis_state(len)], dis_slot_bits );
        if( dis_slot < start_dis_model ) rep0 = dis_slot;
        else
          {
          const int direct_bits = ( dis_slot >> 1 ) - 1;
          rep0 = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
          if( dis_slot < end_dis_model )
            rep0 += range_decoder.decode_tree_reversed( bm_dis + rep0 - dis_slot, direct_bits );
          else
            {
            rep0 += range_decoder.decode( direct_bits - dis_align_bits ) << dis_align_bits;
            rep0 += range_decoder.decode_tree_reversed( bm_align, dis_align_bits );
            if( rep0 == 0xFFFFFFFF )		// Marker found
              {
              rep0 = rep0_saved;
              range_decoder.normalize();
              if( len == min_match_len )	// End Of Stream marker
                {
                member_finished_ = true;
                if( verify_trailer() ) return 0; else return 3;
                }
              if( len == min_match_len + 1 )	// Sync Flush marker
                {
                if( range_decoder.try_reload( true ) ) continue;
                else return 0;
                }
              return 4;
              }
            if( rep0 >= (unsigned int)dictionary_size ) return 1;
            }
          }
        rep3 = rep2; rep2 = rep1; rep1 = rep0_saved;
        state.set_match();
        }
      copy_block( rep0, len );
      }
    }
  }
