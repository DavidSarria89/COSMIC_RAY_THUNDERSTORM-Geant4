// Function to be called from the fortran PARMA code
// to avoid a very nasty bug when reading input file
// where dots (".") are interpreter as separator, or end of number
// and commas (",") as decimal separator
#define _GNU_SOURCE

#include <locale.h>

void
setlocale_to_avoid_bug()
{
    setlocale(LC_ALL, "C");
    return;
}

// void str_to_double(char* sstring, double *dd) {
//  setlocale(LC_ALL, "C");
//  printf ( "%s\n", sstring );
//  *dd = atof(sstring);
//  printf ( "%lf\n", *dd );
// return;
//}
