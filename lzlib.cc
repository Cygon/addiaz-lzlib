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

#include <algorithm>
#include <cstring>
#include <vector>
#include <stdint.h>

#include "lzlib.h"
#include "lzip.h"
#include "decoder.h"
#include "encoder.h"


namespace {

struct Encoder
  {
  long long partial_in_size;
  long long partial_out_size;
  Matchfinder * matchfinder;
  LZ_encoder * lz_encoder;
  LZ_errno lz_errno;
  const File_header member_header;

  Encoder( const File_header & header ) throw()
    :
    partial_in_size( 0 ),
    partial_out_size( 0 ),
    matchfinder( 0 ),
    lz_encoder( 0 ),
    lz_errno( LZ_ok ),
    member_header( header )
    {}
  };


bool verify_encoder( void * const encoder )
  {
  if( !encoder ) return false;
  if( !((Encoder *)encoder)->matchfinder ||
      !((Encoder *)encoder)->lz_encoder )
    { ((Encoder *)encoder)->lz_errno = LZ_bad_argument; return false; }
  return true;
  }


struct Decoder
  {
  long long partial_in_size;
  long long partial_out_size;
  Input_buffer * ibuf;
  LZ_decoder * lz_decoder;
  LZ_errno lz_errno;

  Decoder() throw()
    :
    partial_in_size( 0 ),
    partial_out_size( 0 ),
    ibuf( 0 ),
    lz_decoder( 0 ),
    lz_errno( LZ_ok )
    {}
  };


bool verify_decoder( void * const decoder )
  {
  if( !decoder ) return false;
  if( !((Decoder *)decoder)->ibuf )
    { ((Decoder *)decoder)->lz_errno = LZ_bad_argument; return false; }
  return true;
  }

} // end namespace




const char * LZ_version() { return LZ_version_string; }


void * LZ_compress_open( const int dictionary_size, const int match_len_limit,
                         const long long member_size )
  {
  File_header header;
  header.set_magic();
  const bool error = ( !header.dictionary_size( dictionary_size ) ||
                       match_len_limit < 5 || match_len_limit > max_match_len );

  Encoder * encoder = new( std::nothrow ) Encoder( header );
  if( !encoder ) return 0;
  if( error ) ((Encoder *)encoder)->lz_errno = LZ_bad_argument;
  else
    {
    try {
      encoder->matchfinder =
        new Matchfinder( header.dictionary_size(), match_len_limit );
      }
    catch( std::bad_alloc ) { encoder->matchfinder = 0; }
    if( encoder->matchfinder )
      {
      try {
        encoder->lz_encoder =
          new LZ_encoder( *encoder->matchfinder, header, member_size );
        }
      catch( std::bad_alloc )
        { delete encoder->matchfinder;  encoder->matchfinder = 0;
          encoder->lz_encoder = 0; }
      }
    if( !encoder->lz_encoder ) ((Encoder *)encoder)->lz_errno = LZ_mem_error;
    }
  return encoder;
  }


int LZ_compress_close( void * const encoder )
  {
  if( !encoder ) return -1;
  if( ((Encoder *)encoder)->lz_encoder )
    delete ((Encoder *)encoder)->lz_encoder;
  if( ((Encoder *)encoder)->matchfinder )
    delete ((Encoder *)encoder)->matchfinder;
  delete (Encoder *)encoder;
  return 0;
  }


int LZ_compress_finish( void * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  ((Encoder *)encoder)->matchfinder->finish();
  return 0;
  }


int LZ_compress_finish_member( void * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  ((Encoder *)encoder)->lz_encoder->finish_member();
  return 0;
  }


int LZ_compress_restart_member( void * const encoder,
                                const long long member_size )
  {
  if( !verify_encoder( encoder ) ) return -1;
  if( !((Encoder *)encoder)->lz_encoder->member_finished() )
    { ((Encoder *)encoder)->lz_errno = LZ_sequence_error; return -1; }

  ((Encoder *)encoder)->partial_in_size +=
    ((Encoder *)encoder)->matchfinder->data_position();
  ((Encoder *)encoder)->partial_out_size +=
    ((Encoder *)encoder)->lz_encoder->member_position();
  if( !((Encoder *)encoder)->matchfinder->reset() )
    { ((Encoder *)encoder)->lz_errno = LZ_library_error; return -1; }

  delete ((Encoder *)encoder)->lz_encoder;
  try {
    ((Encoder *)encoder)->lz_encoder =
      new LZ_encoder( *((Encoder *)encoder)->matchfinder,
                      ((Encoder *)encoder)->member_header, member_size );
    }
  catch( std::bad_alloc )
    { ((Encoder *)encoder)->lz_encoder = 0;
      ((Encoder *)encoder)->lz_errno = LZ_mem_error; return -1; }
  return 0;
  }


int LZ_compress_read( void * const encoder, uint8_t * const buffer,
                      const int size )
  {
  if( !verify_encoder( encoder ) ) return -1;
  if( !((Encoder *)encoder)->lz_encoder->encode_member() )
    { ((Encoder *)encoder)->lz_errno = LZ_library_error; return -1; }
  return ((Encoder *)encoder)->lz_encoder->read_data( buffer, size );
  }


int LZ_compress_write( void * const encoder, uint8_t * const buffer,
                       const int size )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return ((Encoder *)encoder)->matchfinder->write_data( buffer, size );
  }


enum LZ_errno LZ_compress_errno( void * const encoder )
  {
  if( !encoder ) return LZ_bad_argument;
  return ((Encoder *)encoder)->lz_errno;
  }


int LZ_compress_finished( void * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return ((Encoder *)encoder)->matchfinder->finished();
  }


int LZ_compress_member_finished( void * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return ((Encoder *)encoder)->lz_encoder->member_finished();
  }


long long LZ_compress_data_position( void * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return ((Encoder *)encoder)->matchfinder->data_position();
  }


long long LZ_compress_member_position( void * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return ((Encoder *)encoder)->lz_encoder->member_position();
  }


long long LZ_compress_total_in_size( void * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return ((Encoder *)encoder)->partial_in_size +
         ((Encoder *)encoder)->matchfinder->data_position();
  }


long long LZ_compress_total_out_size( void * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return ((Encoder *)encoder)->partial_out_size +
         ((Encoder *)encoder)->lz_encoder->member_position();
  }


void * LZ_decompress_open()
  {
  Decoder * decoder = new( std::nothrow ) Decoder();
  if( !decoder ) return 0;

  try { decoder->ibuf = new Input_buffer(); }
  catch( std::bad_alloc )
    { decoder->ibuf = 0; ((Decoder *)decoder)->lz_errno = LZ_mem_error; }
  return decoder;
  }


int LZ_decompress_close( void * const decoder )
  {
  if( !decoder ) return -1;
  if( ((Decoder *)decoder)->lz_decoder )
    delete ((Decoder *)decoder)->lz_decoder;
  if( ((Decoder *)decoder)->ibuf )
    delete ((Decoder *)decoder)->ibuf;
  delete (Decoder *)decoder;
  return 0;
  }


int LZ_decompress_finish( void * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  ((Decoder *)decoder)->ibuf->finish();
  return 0;
  }


int LZ_decompress_read( void * const decoder, uint8_t * const buffer,
                        const int size )
  {
  if( !verify_decoder( decoder ) ) return -1;
  if( !((Decoder *)decoder)->lz_decoder )
    {
    if( ((Decoder *)decoder)->ibuf->used_bytes() <
        5 + (int)sizeof( File_header ) )
      {
      if( ((Decoder *)decoder)->ibuf->at_stream_end() &&
          !((Decoder *)decoder)->ibuf->finished() )	// trailing garbage
        { ((Decoder *)decoder)->lz_errno = LZ_header_error; return -1; }
      else return 0;
      }
    File_header header;
    for( unsigned int i = 0; i < sizeof header; ++i )
      ((uint8_t *)&header)[i] = ((Decoder *)decoder)->ibuf->get_byte();
    if( !header.verify_magic() || !header.verify_version() ||
        header.dictionary_size() < min_dictionary_size ||
        header.dictionary_size() > max_dictionary_size )
      { ((Decoder *)decoder)->lz_errno = LZ_header_error; return -1; }
    try {
      ((Decoder *)decoder)->lz_decoder =
        new LZ_decoder( header, *((Decoder *)decoder)->ibuf );
      }
    catch( std::bad_alloc )
      { ((Decoder *)decoder)->lz_decoder = 0;
        ((Decoder *)decoder)->lz_errno = LZ_mem_error; return -1; }
    }
  const int result = ((Decoder *)decoder)->lz_decoder->decode_member();
  if( result != 0 )
    {
    if( result == 2 ) ((Decoder *)decoder)->lz_errno = LZ_unexpected_eof;
    else ((Decoder *)decoder)->lz_errno = LZ_data_error;
    return -1;
    }
  const int out_size = ((Decoder *)decoder)->lz_decoder->read_data( buffer, size );
  if( ((Decoder *)decoder)->lz_decoder->member_finished() )
    {
    ((Decoder *)decoder)->partial_in_size +=
      ((Decoder *)decoder)->lz_decoder->member_position();
    ((Decoder *)decoder)->partial_out_size +=
      ((Decoder *)decoder)->lz_decoder->data_position();
    delete ((Decoder *)decoder)->lz_decoder;
    ((Decoder *)decoder)->lz_decoder = 0;
    }
  return out_size;
  }


int LZ_decompress_write( void * const decoder, uint8_t * const buffer,
                         const int size )
  {
  if( !verify_decoder( decoder ) ) return -1;
  return ((Decoder *)decoder)->ibuf->write_data( buffer, size );
  }


enum LZ_errno LZ_decompress_errno( void * const decoder )
  {
  if( !decoder ) return LZ_bad_argument;
  return ((Decoder *)decoder)->lz_errno;
  }


int LZ_decompress_finished( void * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  return ((Decoder *)decoder)->ibuf->finished();
  }


long long LZ_decompress_data_position( void * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  if( ((Decoder *)decoder)->lz_decoder )
    return ((Decoder *)decoder)->lz_decoder->data_position();
  else return 0;
  }


long long LZ_decompress_member_position( void * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  if( ((Decoder *)decoder)->lz_decoder )
    return ((Decoder *)decoder)->lz_decoder->member_position();
  else return 0;
  }


long long LZ_decompress_total_in_size( void * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  return ((Decoder *)decoder)->partial_in_size +
         ((Decoder *)decoder)->lz_decoder ?
         ((Decoder *)decoder)->lz_decoder->member_position() : 0;
  }


long long LZ_decompress_total_out_size( void * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  return ((Decoder *)decoder)->partial_out_size +
         ((Decoder *)decoder)->lz_decoder ?
         ((Decoder *)decoder)->lz_decoder->data_position() : 0;
  }
