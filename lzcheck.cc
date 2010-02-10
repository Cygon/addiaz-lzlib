/*  Lzcheck - A test program for the lzlib library
    Copyright (C) 2009, 2010 Antonio Diaz Diaz.

    This program is free software: you have unlimited permission
    to copy, distribute and modify it.

    Usage is:
      lzcheck filename.txt

    This program reads the specified text file and then compresses it,
    line by line, to test the flushing mechanism and the member
    restart/reset/sync functions.
*/

#define _FILE_OFFSET_BITS 64

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <unistd.h>

#include "lzlib.h"

#ifndef LLONG_MAX
#define LLONG_MAX  0x7FFFFFFFFFFFFFFFLL
#endif
#ifndef LLONG_MIN
#define LLONG_MIN  (-LLONG_MAX - 1LL)
#endif
#ifndef ULLONG_MAX
#define ULLONG_MAX 0xFFFFFFFFFFFFFFFFULL
#endif

const int buffer_size = 32768;
uint8_t in_buffer[buffer_size];
uint8_t mid_buffer[buffer_size];
uint8_t out_buffer[buffer_size];


int main( const int argc, const char * argv[] )
  {
  if( argc < 2 )
    {
    std::fprintf( stderr, "Usage: lzcheck filename.txt\n" );
    return 1;
    }

  FILE *file = std::fopen( argv[1], "rb" );
  if( !file )
    {
    std::fprintf( stderr, "Can't open file `%s' for reading\n", argv[1] );
    return 1;
    }
//  std::fprintf( stderr, "lzcheck: testing file `%s'\n", argv[1] );

  const int dictionary_size = 1 << 20;
  const int match_len_limit = 80;
  const long long member_size = LLONG_MAX;
  LZ_Encoder * encoder = LZ_compress_open( dictionary_size, match_len_limit,
                                           member_size );
  if( !encoder || LZ_compress_errno( encoder ) != LZ_ok )
    {
    const bool mem_error = ( LZ_compress_errno( encoder ) == LZ_mem_error );
    LZ_compress_close( encoder );
    if( mem_error )
      {
      std::fprintf( stderr, "not enough memory.\n" );
      return 1;
      }
    std::fprintf( stderr, "internal error: invalid argument to encoder.\n" );
    return 3;
    }

  LZ_Decoder * decoder = LZ_decompress_open();
  if( !decoder || LZ_decompress_errno( decoder ) != LZ_ok )
    {
    LZ_decompress_close( decoder );
    std::fprintf( stderr, "not enough memory.\n" );
    return 1;
    }

  int retval = 0;
  while( retval <= 1 )
    {
    const int read_size = std::fread( in_buffer, 1, buffer_size, file );
    if( read_size <= 0 ) break;			// end of file

    for( int l = 0, r = 1; r <= read_size; l = r, ++r )
      {
      while( r < read_size && in_buffer[r-1] != '\n' ) ++r;
      const int in_size = LZ_compress_write( encoder, in_buffer + l, r - l );
      if( in_size < r - l ) r = l + in_size;
      LZ_compress_sync_flush( encoder );
      const int mid_size = LZ_compress_read( encoder, mid_buffer, buffer_size );
      if( mid_size < 0 )
        {
        std::fprintf( stderr, "LZ_compress_read error: %s.\n",
                      LZ_strerror( LZ_compress_errno( encoder ) ) );
        retval = 3; break;
        }
      LZ_decompress_write( decoder, mid_buffer, mid_size );
      const int out_size = LZ_decompress_read( decoder, out_buffer, buffer_size );
      if( out_size < 0 )
        {
        std::fprintf( stderr, "LZ_decompress_read error: %s.\n",
                      LZ_strerror( LZ_decompress_errno( decoder ) ) );
        retval = 3; break;
        }

      if( out_size != in_size || std::memcmp( in_buffer + l, out_buffer, out_size ) )
        {
        std::fprintf( stderr, "sync error at pos %d. in_size = %d, out_size = %d\n",
                      l, in_size, out_size );
        for( int i = 0; i < in_size; ++i )
          std::fputc( in_buffer[l+i], stderr );
        if( in_buffer[l+in_size-1] != '\n' )
          std::fputc( '\n', stderr );
        for( int i = 0; i < out_size; ++i )
          std::fputc( out_buffer[i], stderr );
        std::fputc( '\n', stderr );
        retval = 1;
        }
      }
    }

  if( retval <= 1 )
    {
    std::rewind( file );
    if( LZ_compress_finish( encoder ) < 0 ||
        LZ_decompress_write( decoder, mid_buffer, LZ_compress_read( encoder, mid_buffer, buffer_size ) ) < 0 ||
        LZ_decompress_read( decoder, out_buffer, buffer_size ) != 0 ||
        LZ_compress_restart_member( encoder, member_size ) < 0 )
      {
      std::fprintf( stderr, "can't finish member: %s.\n",
                    LZ_strerror( LZ_decompress_errno( decoder ) ) );
      retval = 3;
      }
    }

  while( retval <= 1 )
    {
    const int read_size = std::fread( in_buffer, 1, buffer_size / 2, file );
    if( read_size <= 0 ) break;			// end of file

    for( int l = 0, r = 1; r <= read_size; l = r, ++r )
      {
      while( r < read_size && in_buffer[r-1] != '\n' ) ++r;
      const int leading_garbage = (l == 0) ? std::min( r, read_size / 2 ) : 0;
      const int in_size = LZ_compress_write( encoder, in_buffer + l, r - l );
      if( in_size < r - l ) r = l + in_size;
      LZ_compress_sync_flush( encoder );
      if( leading_garbage )
        std::memset( mid_buffer, in_buffer[0], leading_garbage );
      const int mid_size = LZ_compress_read( encoder,
                                             mid_buffer + leading_garbage,
                                             buffer_size - leading_garbage );
      if( mid_size < 0 )
        {
        std::fprintf( stderr, "LZ_compress_read error: %s.\n",
                      LZ_strerror( LZ_compress_errno( encoder ) ) );
        retval = 3; break;
        }
      LZ_decompress_write( decoder, mid_buffer, mid_size + leading_garbage );
      int out_size = LZ_decompress_read( decoder, out_buffer, buffer_size );
      if( out_size < 0 )
        {
        if( LZ_decompress_errno( decoder ) == LZ_header_error ||
            LZ_decompress_errno( decoder ) == LZ_data_error )
          {
          LZ_decompress_sync_to_member( decoder );  // remove leading garbage
          out_size = LZ_decompress_read( decoder, out_buffer, buffer_size );
          }
        if( out_size < 0 )
          {
          std::fprintf( stderr, "LZ_decompress_read error: %s.\n",
                        LZ_strerror( LZ_decompress_errno( decoder ) ) );
          retval = 3; break;
          }
        }

      if( out_size != in_size || std::memcmp( in_buffer + l, out_buffer, out_size ) )
        {
        std::fprintf( stderr, "sync error at pos %d. in_size = %d, out_size = %d\n",
                      l, in_size, out_size );
        for( int i = 0; i < in_size; ++i )
          std::fputc( in_buffer[l+i], stderr );
        if( in_buffer[l+in_size-1] != '\n' )
          std::fputc( '\n', stderr );
        for( int i = 0; i < out_size; ++i )
          std::fputc( out_buffer[i], stderr );
        std::fputc( '\n', stderr );
        retval = 1;
        }
      }
    if( retval >= 3 ) break;

    if( LZ_compress_finish( encoder ) < 0 ||
        LZ_decompress_write( decoder, mid_buffer, LZ_compress_read( encoder, mid_buffer, buffer_size ) ) < 0 ||
        LZ_decompress_read( decoder, out_buffer, buffer_size ) != 0 ||
        LZ_decompress_reset( decoder ) < 0 ||
        LZ_compress_restart_member( encoder, member_size ) < 0 )
      {
      std::fprintf( stderr, "can't restart member: %s.\n",
                    LZ_strerror( LZ_decompress_errno( decoder ) ) );
      retval = 3; break;
      }

    const int size = std::min( 100, read_size );
    if( LZ_compress_write( encoder, in_buffer, size ) != size ||
        LZ_compress_finish( encoder ) < 0 ||
        LZ_decompress_write( decoder, mid_buffer, LZ_compress_read( encoder, mid_buffer, buffer_size ) ) < 0 ||
        LZ_decompress_read( decoder, out_buffer, 0 ) != 0 ||
        LZ_decompress_sync_to_member( decoder ) < 0 ||
        LZ_compress_restart_member( encoder, member_size ) < 0 )
      {
      std::fprintf( stderr, "can't seek to next member: %s.\n",
                    LZ_strerror( LZ_decompress_errno( decoder ) ) );
      retval = 3; break;
      }
    }

  LZ_decompress_close( decoder );
  LZ_compress_close( encoder );
  std::fclose( file );
  return retval;
  }
