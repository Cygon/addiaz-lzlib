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

#ifdef __cplusplus
extern "C" {
#endif

const char * const LZ_version_string = "0.5";

enum { min_dictionary_bits = 12,
       min_dictionary_size = 1 << min_dictionary_bits,
       max_dictionary_bits = 29,
       max_dictionary_size = 1 << max_dictionary_bits };

enum LZ_errno { LZ_ok = 0, LZ_bad_argument, LZ_mem_error, LZ_sequence_error,
                LZ_header_error, LZ_unexpected_eof, LZ_data_error,
                LZ_library_error };


const char * LZ_version( void );


void * LZ_compress_open( const int dictionary_size, const int match_len_limit,
                         const long long member_size );
int LZ_compress_restart_member( void * const encoder,
                                const long long member_size );
int LZ_compress_close( void * const encoder );
int LZ_compress_finish( void * const encoder );
int LZ_compress_sync_flush( void * const encoder );

int LZ_compress_read( void * const encoder, uint8_t * const buffer,
                      const int size );
int LZ_compress_write( void * const encoder, uint8_t * const buffer,
                       const int size );
int LZ_compress_write_size( void * const encoder );

enum LZ_errno LZ_compress_errno( void * const encoder );
int LZ_compress_finished( void * const encoder );
int LZ_compress_member_finished( void * const encoder );

long long LZ_compress_data_position( void * const encoder );
long long LZ_compress_member_position( void * const encoder );
long long LZ_compress_total_in_size( void * const encoder );
long long LZ_compress_total_out_size( void * const encoder );


void * LZ_decompress_open( void );
int LZ_decompress_close( void * const decoder );
int LZ_decompress_finish( void * const decoder );

int LZ_decompress_read( void * const decoder, uint8_t * const buffer,
                        const int size );
int LZ_decompress_write( void * const decoder, uint8_t * const buffer,
                         const int size );

enum LZ_errno LZ_decompress_errno( void * const decoder );
int LZ_decompress_finished( void * const decoder );

long long LZ_decompress_data_position( void * const decoder );
long long LZ_decompress_member_position( void * const decoder );
long long LZ_decompress_total_in_size( void * const decoder );
long long LZ_decompress_total_out_size( void * const decoder );

#ifdef __cplusplus
}
#endif
