/*  Minilzip - A test program for the lzlib library
    Copyright (C) 2009 Antonio Diaz Diaz.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
    Return values: 0 for a normal exit, 1 for environmental problems
    (file not found, invalid flags, I/O errors, etc), 2 to indicate a
    corrupt or invalid input file, 3 for an internal consistency error
    (eg, bug) which caused lzip to panic.
*/

#define _FILE_OFFSET_BITS 64

#include <algorithm>
#include <cerrno>
#include <climits>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fcntl.h>
#include <stdint.h>
#include <unistd.h>
#include <utime.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "arg_parser.h"
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

void show_error( const char * msg, const int errcode = 0, const bool help = false ) throw();
void internal_error( const char * msg );
int readblock( const int fd, char * buf, const int size ) throw();
int writeblock( const int fd, const char * buf, const int size ) throw();


namespace {

const char * invocation_name = 0;
const char * const Program_name    = "Minilzip";
const char * const program_name    = "minilzip";
const char * const program_year    = "2009";

struct { const char * from; const char * to; } const known_extensions[] = {
  { ".lz",  ""     },
  { ".tlz", ".tar" },
  { 0,      0      } };

struct lzma_options
  {
  int dictionary_size;		// 4KiB..512MiB
  int match_len_limit;		// 5..273
  };

enum Mode { m_compress = 0, m_decompress, m_test };

std::string output_filename;
int outhandle = -1;
int verbosity = 0;
bool delete_output_on_interrupt = false;

class Pretty_print
  {
  const char * const stdin_name;
  const unsigned int stdin_name_len;
  unsigned int longest_name;
  std::string name_;
  mutable bool first_post;

public:
  Pretty_print( const std::vector< std::string > & filenames )
    : stdin_name( "(stdin)" ), stdin_name_len( std::strlen( stdin_name ) ),
      longest_name( 0 ), first_post( false )
    {
    for( unsigned int i = 0; i < filenames.size(); ++i )
      {
      const std::string & s = filenames[i];
      const unsigned int len = ( ( s == "-" ) ? stdin_name_len : s.size() );
      if( len > longest_name ) longest_name = len;
      }
    if( longest_name == 0 ) longest_name = stdin_name_len;
    }

  void set_name( const std::string & filename )
    {
    if( filename.size() && filename != "-" ) name_ = filename;
    else name_ = stdin_name;
    first_post = true;
    }

  void reset() const throw() { if( name_.size() ) first_post = true; }
  const char * name() const throw() { return name_.c_str(); }
  void operator()( const char * const msg = 0 ) const throw();
  };


void show_help() throw()
  {
  std::printf( "%s - A test program for the lzlib library.\n", Program_name );
  std::printf( "\nUsage: %s [options] [files]\n", invocation_name );
  std::printf( "\nOptions:\n" );
  std::printf( "  -h, --help                 display this help and exit\n" );
  std::printf( "  -V, --version              output version information and exit\n" );
  std::printf( "  -b, --member-size=<n>      set member size limit in bytes\n" );
  std::printf( "  -c, --stdout               send output to standard output\n" );
  std::printf( "  -d, --decompress           decompress\n" );
  std::printf( "  -f, --force                overwrite existing output files\n" );
  std::printf( "  -k, --keep                 keep (don't delete) input files\n" );
  std::printf( "  -m, --match-length=<n>     set match length limit in bytes [80]\n" );
  std::printf( "  -o, --output=<file>        if reading stdin, place the output into <file>\n" );
  std::printf( "  -q, --quiet                suppress all messages\n" );
  std::printf( "  -s, --dictionary-size=<n>  set dictionary size limit in bytes [8MiB]\n" );
  std::printf( "  -S, --volume-size=<n>      set volume size limit in bytes\n" );
  std::printf( "  -t, --test                 test compressed file integrity\n" );
  std::printf( "  -v, --verbose              be verbose (a 2nd -v gives more)\n" );
  std::printf( "  -1 .. -9                   set compression level [default 6]\n" );
  std::printf( "      --fast                 alias for -1\n" );
  std::printf( "      --best                 alias for -9\n" );
  std::printf( "If no file names are given, %s compresses or decompresses\n", program_name );
  std::printf( "from standard input to standard output.\n" );
  std::printf( "Numbers may be followed by a multiplier: k = kB = 10^3 = 1000,\n" );
  std::printf( "Ki = KiB = 2^10 = 1024, M = 10^6, Mi = 2^20, G = 10^9, Gi = 2^30, etc...\n" );
  std::printf( "\nReport bugs to lzip-bug@nongnu.org\n" );
  std::printf( "Lzip home page: http://www.nongnu.org/lzip/lzip.html\n" );
  }


void show_version() throw()
  {
  std::printf( "%s %s\n", Program_name, PROGVERSION );
  std::printf( "Copyright (C) %s Antonio Diaz Diaz.\n", program_year );
  std::printf( "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n" );
  std::printf( "This is free software: you are free to change and redistribute it.\n" );
  std::printf( "There is NO WARRANTY, to the extent permitted by law.\n" );
  }


const char * format_num( long long num, long long limit = 9999,
                         const int set_prefix = 0 ) throw()
  {
  const char * const si_prefix[8] =
    { "k", "M", "G", "T", "P", "E", "Z", "Y" };
  const char * const binary_prefix[8] =
    { "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi", "Yi" };
  static bool si = false;
  static char buf[16];

  if( set_prefix ) si = ( set_prefix > 0 );
  const int factor = ( si ) ? 1000 : 1024;
  const char * const *prefix = ( si ) ? si_prefix : binary_prefix;
  const char *p = "";
  limit = std::max( 999LL, std::min( 999999LL, limit ) );

  for( int i = 0; i < 8 && ( llabs( num ) > limit ||
       ( llabs( num ) >= factor && num % factor == 0 ) ); ++i )
    { num /= factor; p = prefix[i]; }
  snprintf( buf, sizeof buf, "%lld %s", num, p );
  return buf;
  }


long long getnum( const char * ptr, const int bs,
                  const long long llimit = LLONG_MIN + 1,
                  const long long ulimit = LLONG_MAX ) throw()
  {
  errno = 0;
  char *tail;
  long long result = strtoll( ptr, &tail, 0 );
  if( tail == ptr )
    {
    show_error( "bad or missing numerical argument", 0, true );
    std::exit( 1 );
    }

  if( !errno && tail[0] )
    {
    int factor = ( tail[1] == 'i' ) ? 1024 : 1000;
    int exponent = 0;
    bool bad_multiplier = false;
    switch( tail[0] )
      {
      case ' ': break;
      case 'b': if( bs > 0 ) { factor = bs; exponent = 1; }
                else bad_multiplier = true;
                break;
      case 'Y': exponent = 8; break;
      case 'Z': exponent = 7; break;
      case 'E': exponent = 6; break;
      case 'P': exponent = 5; break;
      case 'T': exponent = 4; break;
      case 'G': exponent = 3; break;
      case 'M': exponent = 2; break;
      case 'K': if( factor == 1024 ) exponent = 1; else bad_multiplier = true;
                break;
      case 'k': if( factor == 1000 ) exponent = 1; else bad_multiplier = true;
                break;
      default: bad_multiplier = true;
      }
    if( bad_multiplier )
      {
      show_error( "bad multiplier in numerical argument", 0, true );
      std::exit( 1 );
      }
    for( int i = 0; i < exponent; ++i )
      {
      if( LLONG_MAX / factor >= llabs( result ) ) result *= factor;
      else { errno = ERANGE; break; }
      }
    }
  if( !errno && ( result < llimit || result > ulimit ) ) errno = ERANGE;
  if( errno )
    {
    show_error( "numerical argument out of limits" );
    std::exit( 1 );
    }
  return result;
  }


int get_dict_size( const char * arg ) throw()
  {
  char *tail;
  int bits = std::strtol( arg, &tail, 0 );
  if( bits >= min_dictionary_bits && bits <= max_dictionary_bits && *tail == 0 )
    return ( 1 << bits );
  return getnum( arg, 0, min_dictionary_size, max_dictionary_size );
  }


int extension_index( const std::string & name ) throw()
  {
  for( int i = 0; known_extensions[i].from; ++i )
    {
    const std::string ext( known_extensions[i].from );
    if( name.size() > ext.size() &&
        name.compare( name.size() - ext.size(), ext.size(), ext ) == 0 )
      return i;
    }
  return -1;
  }


int open_instream( const std::string & name, struct stat * in_statsp,
                   const Mode program_mode, const int eindex,
                   const bool force, const bool to_stdout ) throw()
  {
  int inhandle = -1;
  if( program_mode == m_compress && !force && eindex >= 0 )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "%s: input file `%s' already has `%s' suffix.\n",
                    program_name, name.c_str(),
                    known_extensions[eindex].from );
    }
  else
    {
    inhandle = open( name.c_str(), O_RDONLY );
    if( inhandle < 0 )
      {
      if( verbosity >= 0 )
        std::fprintf( stderr, "%s: Can't open input file `%s': %s.\n",
                      program_name, name.c_str(), std::strerror( errno ) );
      }
    else
      {
      const int i = fstat( inhandle, in_statsp );
      const mode_t & mode = in_statsp->st_mode;
      if( i < 0 || !( S_ISREG( mode ) || ( to_stdout &&
                      ( S_ISFIFO( mode ) || S_ISSOCK( mode ) ||
                        S_ISBLK( mode ) || S_ISCHR( mode ) ) ) ) )
        {
        if( verbosity >= 0 )
          std::fprintf( stderr, "%s: input file `%s' is not a regular file%s.\n",
                        program_name, name.c_str(),
                        to_stdout ? "" : " and `--stdout' was not specified" );
        close( inhandle );
        inhandle = -1;
        }
      }
    }
  return inhandle;
  }


void set_c_outname( const std::string & name, const bool multifile ) throw()
  {
  output_filename = name;
  if( multifile ) output_filename += "00001";
  output_filename += known_extensions[0].from;
  }


void set_d_outname( const std::string & name, const int i ) throw()
  {
  if( i >= 0 )
    {
    const std::string from( known_extensions[i].from );
    if( name.size() > from.size() )
      {
      output_filename.assign( name, 0, name.size() - from.size() );
      output_filename += known_extensions[i].to;
      return;
      }
    }
  output_filename = name; output_filename += ".out";
  if( verbosity >= 0 )
    std::fprintf( stderr, "%s: can't guess original name for `%s' -- using `%s'.\n",
                  program_name, name.c_str(), output_filename.c_str() );
  }


bool open_outstream( const bool force ) throw()
  {
  if( force )
    outhandle = open( output_filename.c_str(), O_CREAT | O_TRUNC | O_WRONLY,
                      S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );
  else outhandle = open( output_filename.c_str(), O_CREAT | O_EXCL | O_WRONLY,
                         S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );
  if( outhandle < 0 )
    {
    if( errno == EEXIST ) outhandle = -2; else outhandle = -1;
    if( verbosity >= 0 )
      {
      if( outhandle == -2 )
        std::fprintf( stderr, "%s: Output file %s already exists, skipping.\n",
                      program_name, output_filename.c_str() );
      else
        std::fprintf( stderr, "%s: Can't create output file `%s': %s.\n",
                      program_name, output_filename.c_str(), std::strerror( errno ) );
      }
    }
  return ( outhandle >= 0 );
  }


bool check_tty( const int inhandle, const Mode program_mode ) throw()
  {
  if( program_mode == m_compress && isatty( outhandle ) )
    {
    show_error( "I won't write compressed data to a terminal.", 0, true );
    return false;
    }
  if( ( program_mode == m_decompress || program_mode == m_test ) &&
      isatty( inhandle ) )
    {
    show_error( "I won't read compressed data from a terminal.", 0, true );
    return false;
    }
  return true;
  }


void cleanup_and_fail( const int retval ) throw()
  {
  if( delete_output_on_interrupt )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "%s: Deleting output file `%s', if it exists.\n",
               program_name, output_filename.c_str() );
    if( outhandle >= 0 ) { close( outhandle ); outhandle = -1; }
    if( std::remove( output_filename.c_str() ) != 0 )
      show_error( "WARNING: deletion of output file (apparently) failed." );
    }
  std::exit( retval );
  }


     // Set permissions, owner and times.
void close_and_set_permissions( const struct stat * in_statsp, int * retvalp )
  {
  int tmp = 0;
  if( in_statsp )
    {
    if( fchmod( outhandle, in_statsp->st_mode ) != 0 ) tmp = 1;
    if( !tmp ) (void)fchown( outhandle, in_statsp->st_uid, in_statsp->st_gid );
    // fchown will in many cases return with EPERM, which can be safely ignored.
    }
  if( close( outhandle ) == 0 ) outhandle = -1;
  else cleanup_and_fail( 1 );
  delete_output_on_interrupt = false;
  if( !in_statsp ) return;
  if( !tmp )
    {
    struct utimbuf t;
    t.actime = in_statsp->st_atime;
    t.modtime = in_statsp->st_mtime;
    tmp = utime( output_filename.c_str(), &t );
    }
  if( tmp )
    {
    if( tmp > *retvalp ) *retvalp = tmp;
    show_error( "I can't change output file attributes." );
    cleanup_and_fail( *retvalp );
    }
  }


bool next_filename()
  {
  const unsigned int len = std::strlen( known_extensions[0].from );
  if( output_filename.size() >= len + 5 )		// "*00001.lz"
    for( int i = output_filename.size() - len - 1, j = 0; j < 5; --i, ++j )
      {
      if( output_filename[i] < '9' ) { ++output_filename[i]; return true; }
      else output_filename[i] = '0';
      }
  return false;
  }


int compress( const long long member_size, const long long volume_size,
              lzma_options encoder_options, const int inhandle,
              const Pretty_print & pp, const struct stat * in_statsp,
              int * retvalp )
  {
  if( verbosity >= 1 ) pp();
  void * encoder = LZ_compress_open( encoder_options.dictionary_size,
                                     encoder_options.match_len_limit,
                                     std::min( member_size, volume_size ) );
  if( !encoder || LZ_compress_errno( encoder ) != LZ_ok )
    {
    const bool mem_error = ( LZ_compress_errno( encoder ) == LZ_mem_error );
    LZ_compress_close( encoder );
    if( mem_error )
      { pp( "not enough memory. Try a smaller dictionary size" ); return 1; }
    internal_error( "invalid argument to encoder" );
    }

  long long partial_volume_size = 0;
  const int out_buffer_size = 65536, in_buffer_size = 8 * out_buffer_size;
  uint8_t in_buffer[in_buffer_size], out_buffer[out_buffer_size];
  while( true )
    {
    int in_size = std::min( LZ_compress_write_size( encoder ), in_buffer_size );
    if( in_size > 0 )
      {
      in_size = readblock( inhandle, (char *)in_buffer, in_size );
      if( in_size == 0 ) LZ_compress_finish( encoder );
      else if( in_size != LZ_compress_write( encoder, in_buffer, in_size ) )
        internal_error( "library error" );
//      for( int i = 0; i < 10000; ++i )
//        LZ_compress_sync_flush( encoder );
      }
    int out_size = LZ_compress_read( encoder, out_buffer, out_buffer_size );
//    std::fprintf( stderr, "%6d in_size, %5d out_size.\n", in_size, out_size );
    if( out_size < 0 )
      { pp(); show_error( "read error", errno ); return 1; }
    else if( out_size > 0 )
      {
      const int wr = writeblock( outhandle, (char *)out_buffer, out_size );
      if( wr != out_size )
        { pp(); show_error( "write error", errno ); return 1; }
      }
    else if( in_size == 0 ) internal_error( "library error" );
    if( LZ_compress_member_finished( encoder ) )
      {
      if( LZ_compress_finished( encoder ) == 1 ) break;
      partial_volume_size += LZ_compress_member_position( encoder );
      if( partial_volume_size >= volume_size - min_dictionary_size )
        {
        partial_volume_size = 0;
        if( delete_output_on_interrupt )
          {
          close_and_set_permissions( in_statsp, retvalp );
          if( !next_filename() )
            { pp(); show_error( "too many volume files" ); return 1; }
          if( !open_outstream( true ) ) return 1;
          delete_output_on_interrupt = true;
          }
        }
      const long long size =
        std::min( member_size, volume_size - partial_volume_size );
      if( LZ_compress_restart_member( encoder, size ) < 0 )
        { pp(); show_error( "read error", errno ); return 1; }
      }
    }

  if( verbosity >= 1 )
    {
    const long long in_size = LZ_compress_total_in_size( encoder );
    const long long out_size = LZ_compress_total_out_size( encoder );
    if( in_size <= 0 || out_size <= 0 )
      std::fprintf( stderr, "no data compressed.\n" );
    else
      std::fprintf( stderr, "%6.3f:1, %6.3f bits/byte, "
                            "%5.2f%% saved, %lld in, %lld out.\n",
                    (double)in_size / out_size,
                    ( 8.0 * out_size ) / in_size,
                    100.0 * ( 1.0 - ( (double)out_size / in_size ) ),
                    in_size, out_size );
    }
  LZ_compress_close( encoder );
  return 0;
  }


int decompress( const int inhandle, const Pretty_print & pp,
                const bool testing )
  {
  void * decoder = LZ_decompress_open();
  if( !decoder || LZ_decompress_errno( decoder ) != LZ_ok )
    {
    LZ_decompress_close( decoder );
    pp( "not enough memory. Find a machine with more memory" );
    return 1;
    }
  if( verbosity >= 1 ) pp();

  const int in_buffer_size = 65536, out_buffer_size = 8 * in_buffer_size;
  uint8_t in_buffer[in_buffer_size], out_buffer[out_buffer_size];
  int in_pos = 0, in_stream_pos = 0;
  bool finished = false;
  while( true )
    {
    int in_size = 0;
    if( !finished )
      {
      if( in_stream_pos == 0 )
        in_stream_pos = readblock( inhandle, (char *)in_buffer, in_buffer_size );
      if( in_pos < in_stream_pos )
        {
        in_size = LZ_decompress_write( decoder, in_buffer + in_pos, in_stream_pos - in_pos );
        in_pos += in_size;
        }
      if( in_pos >= in_stream_pos )
        {
        if( in_stream_pos < in_buffer_size )
          { finished = true; LZ_decompress_finish( decoder ); }
        in_stream_pos = 0; in_pos = 0;
        }
      }
    int out_size = LZ_decompress_read( decoder, out_buffer, out_buffer_size );
//    std::fprintf( stderr, "%5d in_size, %6d out_size.\n", in_size, out_size );
    if( out_size < 0 )
      {
      const LZ_errno lz_errno = LZ_decompress_errno( decoder );
      if( lz_errno == LZ_header_error )
        {
        if( LZ_decompress_total_out_size( decoder ) > 0 )
          break;				// trailing garbage
        pp( "error reading member header" );
        return 1;
        }
      if( lz_errno == LZ_mem_error )
        {
        pp( "not enough memory. Find a machine with more memory" );
        return 1;
        }
      if( lz_errno == LZ_unexpected_eof )
        {
        if( verbosity >= 0 )
          { pp();
            std::fprintf( stderr, "file ends unexpectedly at pos %lld\n",
                          LZ_decompress_total_in_size( decoder ) ); }
        return 2;
        }
      pp(); show_error( "read error", errno ); return 1;
      }
    else if( out_size > 0 && outhandle >= 0 )
      {
      const int wr = writeblock( outhandle, (char *)out_buffer, out_size );
      if( wr != out_size )
        { pp(); show_error( "write error", errno ); return 1; }
      }
    if( LZ_decompress_finished( decoder ) == 1 ) break;
    if( finished && in_size == 0 && out_size == 0 )
      internal_error( "library error" );
    }
  if( verbosity >= 1 )
    { if( testing ) std::fprintf( stderr, "ok\n" );
      else std::fprintf( stderr, "done\n" ); }
  LZ_decompress_close( decoder );
  return 0;
  }


extern "C" void signal_handler( int ) throw()
  {
  show_error( "Control-C or similar caught, quitting." );
  cleanup_and_fail( 0 );
  }


void set_signals() throw()
  {
  signal( SIGTERM, signal_handler );
  signal( SIGHUP, signal_handler );
  signal( SIGINT, signal_handler );
  }

} // end namespace


void Pretty_print::operator()( const char * const msg ) const throw()
  {
  if( verbosity >= 0 )
    {
    if( first_post )
      {
      first_post = false;
      std::fprintf( stderr, "  %s: ", name_.c_str() );
      for( unsigned int i = 0; i < longest_name - name_.size(); ++i )
        std::fprintf( stderr, " " );
      if( !msg ) std::fflush( stderr );
      }
    if( msg ) std::fprintf( stderr, "%s.\n", msg );
    }
  }


void show_error( const char * msg, const int errcode, const bool help ) throw()
  {
  if( verbosity >= 0 )
    {
    if( msg && msg[0] != 0 )
      {
      std::fprintf( stderr, "%s: %s", program_name, msg );
      if( errcode > 0 ) std::fprintf( stderr, ": %s", std::strerror( errcode ) );
      std::fprintf( stderr, "\n" );
      }
    if( help && invocation_name && invocation_name[0] != 0 )
      std::fprintf( stderr, "Try `%s --help' for more information.\n", invocation_name );
    }
  }


void internal_error( const char * msg )
  {
  std::string s( "internal error: " ); s += msg;
  show_error( s.c_str() );
  std::exit( 3 );
  }


// Returns the number of bytes really read.
// If (returned value < size) and (errno == 0), means EOF was reached.
//
int readblock( const int fd, char * buf, const int size ) throw()
  {
  int rest = size;
  errno = 0;
  while( rest > 0 )
    {
    errno = 0;
    const int n = read( fd, buf + size - rest, rest );
    if( n > 0 ) rest -= n;
    else if( n == 0 ) break;
    else if( errno != EINTR && errno != EAGAIN ) break;
    }
  return ( rest > 0 ) ? size - rest : size;
  }


// Returns the number of bytes really written.
// If (returned value < size), it is always an error.
//
int writeblock( const int fd, const char * buf, const int size ) throw()
  {
  int rest = size;
  errno = 0;
  while( rest > 0 )
    {
    errno = 0;
    const int n = write( fd, buf + size - rest, rest );
    if( n > 0 ) rest -= n;
    else if( errno && errno != EINTR && errno != EAGAIN ) break;
    }
  return ( rest > 0 ) ? size - rest : size;
  }


int main( const int argc, const char * argv[] )
  {
  // Mapping from gzip/bzip2 style 1..9 compression modes
  // to the corresponding LZMA compression modes.
  const lzma_options option_mapping[] =
    {
    { 1 << 20,  10 },		// -1
    { 1 << 20,  12 },		// -2
    { 1 << 20,  17 },		// -3
    { 1 << 21,  26 },		// -4
    { 1 << 22,  44 },		// -5
    { 1 << 23,  80 },		// -6
    { 1 << 24, 108 },		// -7
    { 1 << 24, 163 },		// -8
    { 1 << 25, 273 } };		// -9
  lzma_options encoder_options = option_mapping[5];	// default = "-6"
  long long member_size = LLONG_MAX;
  long long volume_size = LLONG_MAX;
  int inhandle = -1;
  Mode program_mode = m_compress;
  bool force = false;
  bool keep_input_files = false;
  bool to_stdout = false;
  std::string input_filename;
  std::string default_output_filename;
  std::vector< std::string > filenames;
  invocation_name = argv[0];

  if( LZ_version()[0] != LZ_version_string[0] )
    internal_error( "bad library version" );

  if( std::strcmp( PROGVERSION, LZ_version_string ) )
    internal_error( "bad library version_string" );

  const Arg_parser::Option options[] =
    {
    { '1', "fast",            Arg_parser::no  },
    { '2',  0,                Arg_parser::no  },
    { '3',  0,                Arg_parser::no  },
    { '4',  0,                Arg_parser::no  },
    { '5',  0,                Arg_parser::no  },
    { '6',  0,                Arg_parser::no  },
    { '7',  0,                Arg_parser::no  },
    { '8',  0,                Arg_parser::no  },
    { '9', "best",            Arg_parser::no  },
    { 'b', "member-size",     Arg_parser::yes },
    { 'c', "stdout",          Arg_parser::no  },
    { 'd', "decompress",      Arg_parser::no  },
    { 'f', "force",           Arg_parser::no  },
    { 'h', "help",            Arg_parser::no  },
    { 'k', "keep",            Arg_parser::no  },
    { 'm', "match-length",    Arg_parser::yes },
    { 'o', "output",          Arg_parser::yes },
    { 'q', "quiet",           Arg_parser::no  },
    { 's', "dictionary-size", Arg_parser::yes },
    { 'S', "volume-size",     Arg_parser::yes },
    { 't', "test",            Arg_parser::no  },
    { 'v', "verbose",         Arg_parser::no  },
    { 'V', "version",         Arg_parser::no  },
    {  0 ,  0,                Arg_parser::no  } };

  Arg_parser parser( argc, argv, options );
  if( parser.error().size() )				// bad option
    { show_error( parser.error().c_str(), 0, true ); return 1; }

  int argind = 0;
  for( ; argind < parser.arguments(); ++argind )
    {
    const int code = parser.code( argind );
    if( !code ) break;					// no more options
    const char * arg = parser.argument( argind ).c_str();
    switch( code )
      {
      case '1': case '2': case '3':
      case '4': case '5': case '6':
      case '7': case '8': case '9':
                encoder_options = option_mapping[code-'1']; break;
      case 'b': member_size = getnum( arg, 0, 100000, LLONG_MAX / 2 ); break;
      case 'c': to_stdout = true; break;
      case 'd': program_mode = m_decompress; break;
      case 'f': force = true; break;
      case 'h': show_help(); return 0;
      case 'k': keep_input_files = true; break;
      case 'm': encoder_options.match_len_limit =
                getnum( arg, 0, 5, 273 ); break;
      case 'o': default_output_filename = arg; break;
      case 'q': verbosity = -1; break;
      case 's': encoder_options.dictionary_size = get_dict_size( arg );
                break;
      case 'S': volume_size = getnum( arg, 0, 100000, LLONG_MAX / 2 ); break;
      case 't': program_mode = m_test; break;
      case 'v': if( verbosity < 4 ) ++verbosity; break;
      case 'V': show_version(); return 0;
      default : internal_error( "uncaught option" );
      }
    }

  bool filenames_given = false;
  for( ; argind < parser.arguments(); ++argind )
    {
    if( parser.argument( argind ) != "-" ) filenames_given = true;
    filenames.push_back( parser.argument( argind ) );
    }

  if( filenames.empty() ) filenames.push_back("-");
  if( filenames_given ) set_signals();

  Pretty_print pp( filenames );
  if( program_mode == m_test )
    outhandle = -1;

  int retval = 0;
  for( unsigned int i = 0; i < filenames.size(); ++i )
    {
    struct stat in_stats;
    output_filename.clear();

    if( !filenames[i].size() || filenames[i] == "-" )
      {
      input_filename.clear();
      inhandle = STDIN_FILENO;
      if( program_mode != m_test )
        {
        if( to_stdout || !default_output_filename.size() )
          outhandle = STDOUT_FILENO;
        else
          {
          if( program_mode == m_compress )
            set_c_outname( default_output_filename, volume_size != LLONG_MAX );
          else output_filename = default_output_filename;
          if( !open_outstream( force ) )
            {
            if( outhandle == -1 && retval < 1 ) retval = 1;
            close( inhandle ); inhandle = -1;
            continue;
            }
          }
        }
      }
    else
      {
      input_filename = filenames[i];
      const int eindex = extension_index( input_filename );
      inhandle = open_instream( input_filename, &in_stats, program_mode,
                                eindex, force, to_stdout );
      if( inhandle < 0 ) continue;
      if( program_mode != m_test )
        {
        if( to_stdout ) outhandle = STDOUT_FILENO;
        else
          {
          if( program_mode == m_compress )
            set_c_outname( input_filename, volume_size != LLONG_MAX );
          else set_d_outname( input_filename, eindex );
          if( !open_outstream( force ) )
            {
            if( outhandle == -1 && retval < 1 ) retval = 1;
            close( inhandle ); inhandle = -1;
            continue;
            }
          }
        }
      }

    if( !check_tty( inhandle, program_mode ) ) return 1;

    if( output_filename.size() && !to_stdout && program_mode != m_test )
      delete_output_on_interrupt = true;
    const struct stat * in_statsp = input_filename.size() ? &in_stats : 0;
    pp.set_name( input_filename );
    int tmp = 0;
    if( program_mode == m_compress )
      tmp = compress( member_size, volume_size, encoder_options, inhandle,
                      pp, in_statsp, &retval );
    else
      tmp = decompress( inhandle, pp, program_mode == m_test );
    if( tmp > retval ) retval = tmp;
    if( tmp && program_mode != m_test ) cleanup_and_fail( retval );

    if( delete_output_on_interrupt )
      close_and_set_permissions( in_statsp, &retval );
    if( input_filename.size() )
      {
      close( inhandle ); inhandle = -1;
      if( !keep_input_files && !to_stdout && program_mode != m_test )
        std::remove( input_filename.c_str() );
      }
    }
  if( outhandle >= 0 ) close( outhandle );
  return retval;
  }
