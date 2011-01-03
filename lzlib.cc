/*  Lzlib - A compression library for lzip files
    Copyright (C) 2009, 2010, 2011 Antonio Diaz Diaz.

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

#include <algorithm>
#include <cstring>
#include <stdint.h>

#include "lzlib.h"
#include "lzip.h"
#include "decoder.h"
#include "encoder.h"


using namespace Lzlib;

struct LZ_Encoder
  {
  long long partial_in_size;
  long long partial_out_size;
  Matchfinder * matchfinder;
  LZ_encoder * lz_encoder;
  LZ_Errno lz_errno;
  int flush_pending;
  const File_header member_header;
  bool fatal;

  LZ_Encoder( const File_header & header ) throw()
    :
    partial_in_size( 0 ),
    partial_out_size( 0 ),
    matchfinder( 0 ),
    lz_encoder( 0 ),
    lz_errno( LZ_ok ),
    flush_pending( 0 ),
    member_header( header ),
    fatal( false )
    {}
  };


struct LZ_Decoder
  {
  long long partial_in_size;
  long long partial_out_size;
  Range_decoder * rdec;
  LZ_decoder * lz_decoder;
  LZ_Errno lz_errno;
  File_header member_header;		// header of current member
  bool fatal;
  bool seeking;

  LZ_Decoder() throw()
    :
    partial_in_size( 0 ),
    partial_out_size( 0 ),
    rdec( 0 ),
    lz_decoder( 0 ),
    lz_errno( LZ_ok ),
    fatal( false ),
    seeking( false )
    {
    for( int i = 0; i < File_header::size; ++i ) member_header.data[i] = 0;
    }
  };


namespace Lzlib {

bool verify_encoder( LZ_Encoder * const encoder )
  {
  if( !encoder ) return false;
  if( !encoder->matchfinder || !encoder->lz_encoder )
    { encoder->lz_errno = LZ_bad_argument; return false; }
  return true;
  }


bool verify_decoder( struct LZ_Decoder * const decoder )
  {
  if( !decoder ) return false;
  if( !decoder->rdec )
    { decoder->lz_errno = LZ_bad_argument; return false; }
  return true;
  }

} // end namespace Lzlib


const char * LZ_version() { return LZ_version_string; }


const char * LZ_strerror( const LZ_Errno lz_errno )
  {
  switch( lz_errno )
    {
    case LZ_ok            : return "ok";
    case LZ_bad_argument  : return "bad argument";
    case LZ_mem_error     : return "not enough memory";
    case LZ_sequence_error: return "sequence error";
    case LZ_header_error  : return "header error";
    case LZ_unexpected_eof: return "unexpected eof";
    case LZ_data_error    : return "data error";
    case LZ_library_error : return "library error";
    }
  return "invalid error code";
  }


int LZ_min_dictionary_bits() { return min_dictionary_bits; }
int LZ_min_dictionary_size() { return min_dictionary_size; }
int LZ_max_dictionary_bits() { return max_dictionary_bits; }
int LZ_max_dictionary_size() { return max_dictionary_size; }
int LZ_min_match_len_limit() { return min_match_len_limit; }
int LZ_max_match_len_limit() { return max_match_len; }


/*---------------------- Compression Functions ----------------------*/

LZ_Encoder * LZ_compress_open( const int dictionary_size,
                               const int match_len_limit,
                               const long long member_size )
  {
  File_header header;
  header.set_magic();
  const bool error = ( !header.dictionary_size( dictionary_size ) ||
                       match_len_limit < min_match_len_limit ||
                       match_len_limit > max_match_len );

  LZ_Encoder * encoder = new( std::nothrow ) LZ_Encoder( header );
  if( !encoder ) return 0;
  LZ_Encoder & e = *encoder;
  if( error ) e.lz_errno = LZ_bad_argument;
  else
    {
    try {
      e.matchfinder = new Matchfinder( header.dictionary_size(), match_len_limit );
      }
    catch( std::bad_alloc ) { e.matchfinder = 0; }
    if( e.matchfinder )
      {
      try {
        e.lz_encoder = new LZ_encoder( *e.matchfinder, header, member_size );
        }
      catch( std::bad_alloc )
        {
        delete e.matchfinder;
        e.matchfinder = 0;
        e.lz_encoder = 0;
        }
      }
    if( !e.lz_encoder ) e.lz_errno = LZ_mem_error;
    }
  if( e.lz_errno != LZ_ok ) e.fatal = true;
  return encoder;
  }


int LZ_compress_close( LZ_Encoder * const encoder )
  {
  if( !encoder ) return -1;
  if( encoder->lz_encoder ) delete encoder->lz_encoder;
  if( encoder->matchfinder ) delete encoder->matchfinder;
  delete encoder;
  return 0;
  }


int LZ_compress_finish( LZ_Encoder * const encoder )
  {
  if( !verify_encoder( encoder ) || encoder->fatal ) return -1;
  encoder->matchfinder->flushing( true );
  encoder->flush_pending = 0;
  return 0;
  }


int LZ_compress_restart_member( LZ_Encoder * const encoder,
                                const long long member_size )
  {
  if( !verify_encoder( encoder ) || encoder->fatal ) return -1;
  LZ_Encoder & e = *encoder;
  if( !e.lz_encoder->member_finished() )
    { e.lz_errno = LZ_sequence_error; return -1; }

  e.partial_in_size += e.matchfinder->data_position();
  e.partial_out_size += e.lz_encoder->member_position();
  e.matchfinder->reset();

  delete e.lz_encoder;
  try {
    e.lz_encoder = new LZ_encoder( *e.matchfinder, e.member_header, member_size );
    }
  catch( std::bad_alloc )
    { e.lz_encoder = 0; e.lz_errno = LZ_mem_error; e.fatal = true; return -1; }
  e.lz_errno = LZ_ok;
  return 0;
  }


int LZ_compress_sync_flush( LZ_Encoder * const encoder )
  {
  if( !verify_encoder( encoder ) || encoder->fatal ) return -1;
  LZ_Encoder & e = *encoder;
  if( !e.flush_pending && !e.matchfinder->at_stream_end() )
    {
    e.flush_pending = 2;	// 2 consecutive markers guarantee decoding
    e.matchfinder->flushing( true );
    if( !e.lz_encoder->encode_member( false ) )
      { e.lz_errno = LZ_library_error; e.fatal = true; return -1; }
    while( e.flush_pending > 0 && e.lz_encoder->sync_flush() )
      { if( --e.flush_pending <= 0 ) e.matchfinder->flushing( false ); }
    }
  return 0;
  }


int LZ_compress_read( LZ_Encoder * const encoder,
                      uint8_t * const buffer, const int size )
  {
  if( !verify_encoder( encoder ) || encoder->fatal ) return -1;
  LZ_Encoder & e = *encoder;
  if( !e.lz_encoder->encode_member( !e.flush_pending ) )
    { e.lz_errno = LZ_library_error; e.fatal = true; return -1; }
  while( e.flush_pending > 0 && e.lz_encoder->sync_flush() )
    { if( --e.flush_pending <= 0 ) e.matchfinder->flushing( false ); }
  return e.lz_encoder->read_data( buffer, size );
  }


int LZ_compress_write( LZ_Encoder * const encoder,
                       const uint8_t * const buffer, const int size )
  {
  if( !verify_encoder( encoder ) || encoder->fatal ) return -1;
  if( encoder->flush_pending ) return 0;
  return encoder->matchfinder->write_data( buffer, size );
  }


int LZ_compress_write_size( LZ_Encoder * const encoder )
  {
  if( !verify_encoder( encoder ) || encoder->fatal ) return -1;
  if( encoder->flush_pending ) return 0;
  return encoder->matchfinder->free_bytes();
  }


LZ_Errno LZ_compress_errno( LZ_Encoder * const encoder )
  {
  if( !encoder ) return LZ_bad_argument;
  return encoder->lz_errno;
  }


int LZ_compress_finished( LZ_Encoder * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return ( !encoder->flush_pending && encoder->matchfinder->finished() &&
           encoder->lz_encoder->member_finished() );
  }


int LZ_compress_member_finished( LZ_Encoder * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return encoder->lz_encoder->member_finished();
  }


long long LZ_compress_data_position( LZ_Encoder * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return encoder->matchfinder->data_position();
  }


long long LZ_compress_member_position( LZ_Encoder * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return encoder->lz_encoder->member_position();
  }


long long LZ_compress_total_in_size( LZ_Encoder * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return encoder->partial_in_size + encoder->matchfinder->data_position();
  }


long long LZ_compress_total_out_size( LZ_Encoder * const encoder )
  {
  if( !verify_encoder( encoder ) ) return -1;
  return encoder->partial_out_size + encoder->lz_encoder->member_position();
  }


/*--------------------- Decompression Functions ---------------------*/

struct LZ_Decoder * LZ_decompress_open()
  {
  LZ_Decoder * decoder = new( std::nothrow ) LZ_Decoder;
  if( !decoder ) return 0;

  LZ_Decoder & d = *decoder;
  try { d.rdec = new Range_decoder; }
  catch( std::bad_alloc )
    { d.rdec = 0; d.lz_errno = LZ_mem_error; d.fatal = true; }
  return decoder;
  }


int LZ_decompress_close( struct LZ_Decoder * const decoder )
  {
  if( !decoder ) return -1;
  if( decoder->lz_decoder ) delete decoder->lz_decoder;
  if( decoder->rdec ) delete decoder->rdec;
  delete decoder;
  return 0;
  }


int LZ_decompress_finish( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) || decoder->fatal ) return -1;
  LZ_Decoder & d = *decoder;
  if( d.seeking ) { d.seeking = false; d.rdec->purge(); }
  else d.rdec->finish();
  return 0;
  }


int LZ_decompress_reset( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  LZ_Decoder & d = *decoder;
  if( d.lz_decoder ) { delete d.lz_decoder; d.lz_decoder = 0; }
  d.partial_in_size = 0;
  d.partial_out_size = 0;
  d.rdec->reset();
  d.lz_errno = LZ_ok;
  d.fatal = false;
  d.seeking = false;
  return 0;
  }


int LZ_decompress_sync_to_member( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  LZ_Decoder & d = *decoder;
  if( d.lz_decoder ) { delete d.lz_decoder; d.lz_decoder = 0; }
  if( d.rdec->find_header() ) d.seeking = false;
  else
    {
    if( !d.rdec->at_stream_end() ) d.seeking = true;
    else { d.seeking = false; d.rdec->purge(); }
    }
  d.lz_errno = LZ_ok;
  d.fatal = false;
  return 0;
  }


int LZ_decompress_read( struct LZ_Decoder * const decoder,
                        uint8_t * const buffer, const int size )
  {
  if( !verify_decoder( decoder ) || decoder->fatal ) return -1;
  LZ_Decoder & d = *decoder;
  if( d.seeking ) return 0;
  if( d.lz_decoder && d.lz_decoder->member_finished() )
    {
    d.partial_in_size += d.rdec->member_position();
    d.partial_out_size += d.lz_decoder->data_position();
    delete d.lz_decoder;
    d.lz_decoder = 0;
    }
  if( !d.lz_decoder )
    {
    if( d.rdec->used_bytes() < 5 + File_header::size )
      {
      if( !d.rdec->at_stream_end() || d.rdec->finished() ) return 0;
      d.rdec->purge();			// remove trailing garbage
      d.lz_errno = LZ_header_error;
      d.fatal = true;
      return -1;
      }
    if( !d.rdec->read_header( d.member_header ) )
      {
      d.lz_errno = LZ_header_error;
      d.fatal = true;
      return -1;
      }
    try { d.lz_decoder = new LZ_decoder( d.member_header, *d.rdec ); }
    catch( std::bad_alloc )		// not enough free memory
      {
      d.lz_decoder = 0;
      d.lz_errno = LZ_mem_error;
      d.fatal = true;
      return -1;
      }
    }
  const int result = d.lz_decoder->decode_member();
  if( result != 0 )
    {
    if( result == 2 ) d.lz_errno = LZ_unexpected_eof;
    else d.lz_errno = LZ_data_error;
    d.fatal = true;
    return -1;
    }
  return d.lz_decoder->read_data( buffer, size );
  }


int LZ_decompress_write( struct LZ_Decoder * const decoder,
                         const uint8_t * const buffer, const int size )
  {
  if( !verify_decoder( decoder ) || decoder->fatal ) return -1;
  LZ_Decoder & d = *decoder;
  int result = d.rdec->write_data( buffer, size );
  while( d.seeking )
    {
    if( d.rdec->find_header() ) d.seeking = false;
    if( result >= size ) break;
    const int size2 = d.rdec->write_data( buffer + result, size - result );
    if( size2 > 0 ) result += size2;
    else break;
    }
  return result;
  }


int LZ_decompress_write_size( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) || decoder->fatal ) return -1;
  return decoder->rdec->free_bytes();
  }


LZ_Errno LZ_decompress_errno( struct LZ_Decoder * const decoder )
  {
  if( !decoder ) return LZ_bad_argument;
  return decoder->lz_errno;
  }


int LZ_decompress_finished( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  return ( decoder->rdec->finished() &&
           ( !decoder->lz_decoder || decoder->lz_decoder->member_finished() ) );
  }


int LZ_decompress_member_finished( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  return ( decoder->lz_decoder && decoder->lz_decoder->member_finished() );
  }


int LZ_decompress_member_version( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  return decoder->member_header.version();
  }


int LZ_decompress_dictionary_size( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  return decoder->member_header.dictionary_size();
  }


unsigned int LZ_decompress_data_crc( struct LZ_Decoder * const decoder )
  {
  if( verify_decoder( decoder ) && decoder->lz_decoder )
    return decoder->lz_decoder->crc();
  else return 0;
  }


long long LZ_decompress_data_position( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  if( decoder->lz_decoder )
    return decoder->lz_decoder->data_position();
  else return 0;
  }


long long LZ_decompress_member_position( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  if( decoder->lz_decoder )
    return decoder->rdec->member_position();
  else return 0;
  }


long long LZ_decompress_total_in_size( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  if( decoder->lz_decoder )
    return decoder->partial_in_size + decoder->rdec->member_position();
  return decoder->partial_in_size;
  }


long long LZ_decompress_total_out_size( struct LZ_Decoder * const decoder )
  {
  if( !verify_decoder( decoder ) ) return -1;
  if( decoder->lz_decoder )
    return decoder->partial_out_size + decoder->lz_decoder->data_position();
  return decoder->partial_out_size;
  }
