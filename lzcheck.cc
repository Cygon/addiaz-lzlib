/*  Lzcheck - A test program for the lzlib library
    Copyright (C) 2009 Antonio Diaz Diaz.

    This program is free software: you have unlimited permission
    to copy, distribute and modify it.

    Usage is:
      lzcheck filename.txt

    This program reads the specified text file and then compresses it,
    line by line, to test the flushing mechanism.
*/

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


int main( const int argc, const char * argv[] )
  {
  if( argc < 2 )
    {
    std::fprintf( stderr, "Usage: lzcheck filename.txt\n" );
    return 1;
    }

  FILE *f = std::fopen( argv[1], "rb" );
  if( !f )
    {
    std::fprintf( stderr, "Can't open file `%s' for reading\n", argv[1] );
    return 1;
    }

  const int in_buffer_size = 1 << 20;
  const int mid_buffer_size = 65536;
  const int out_buffer_size = in_buffer_size;
  uint8_t in_buffer[in_buffer_size];
  uint8_t mid_buffer[mid_buffer_size];
  uint8_t out_buffer[out_buffer_size];
  const int in_size = std::fread( in_buffer, 1, in_buffer_size, f );
  if( in_size >= in_buffer_size )
    {
    std::fprintf( stderr, "input file `%s' too big.\n", argv[1] );
    return 1;
    }
  std::fclose( f );

  const int dictionary_size = in_buffer_size;
  const int match_len_limit = 80;
  const long long member_size = LLONG_MAX;
  void * encoder = LZ_compress_open( dictionary_size, match_len_limit,
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

  void * decoder = LZ_decompress_open();
  if( !decoder || LZ_decompress_errno( decoder ) != LZ_ok )
    {
    LZ_decompress_close( decoder );
    std::fprintf( stderr, "not enough memory.\n" );
    return 1;
    }

  for( int l = 0, r = 0; r < in_size; l = r )
    {
    while( ++r < in_size && in_buffer[r-1] != '\n' ) ;
    LZ_compress_write( encoder, in_buffer + l, r - l );
    LZ_compress_sync_flush( encoder );
    int mid_size = LZ_compress_read( encoder, mid_buffer, mid_buffer_size );
    LZ_decompress_write( decoder, mid_buffer, mid_size );
    int out_size = LZ_decompress_read( decoder, out_buffer, out_buffer_size );

    if( out_size != r - l || std::memcmp( in_buffer + l, out_buffer, out_size ) )
      {
      std::printf( "sync error at pos %d. in_size = %d, out_size = %d\n",
                   l, r - l, out_size );
      for( int i = l; i < r; ++i ) std::putchar( in_buffer[i] );
      if( in_buffer[r-1] != '\n' ) std::putchar( '\n' );
      for( int i = 0; i < out_size; ++i ) std::putchar( out_buffer[i] );
      std::putchar( '\n' );
      }
    }

  LZ_decompress_close( decoder );
  LZ_compress_close( encoder );
  return 0;
  }
