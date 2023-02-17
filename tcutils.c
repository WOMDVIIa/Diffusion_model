#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include "tcutils.h"

/************************************************
*   Get path from the values of
*   Thermo-Calc environment variables
*   TC3_HOME (Thermo-Calc 3.0)
*   TCPATH Older versions and fallback
************************************************/

void getThermoCalcEnvironmentPath(char* pathBuffer)
{
  char* path;
  const char *name1=TC3_HOME;
  const char *name2=TCPATH;
  path = getenv(name1);

    if (path == NULL )
    {
        path = getenv(name2);
	if (path == NULL )
	  {
	    strcpy(pathBuffer, " ");
	  }
    }
    strcpy(pathBuffer, path);
}
//------------------------------------------------------------------------------

/************************************************
*   Get path to temp directory
*   If it can't find it - default
*   to current working directory
************************************************/
void getTempEnvironmentPath(char* pathBuffer)
{
  const char *name1=TEMP;
  char* path;

#ifdef WIN32
    path = getenv(name1);
#else
    path=name1;
    strcpy(pathBuffer,name1);
#endif
    if (path == NULL )
    {
        getCurrentWorkingDir(pathBuffer, FILENAME_MAX);
        return;
    }
#ifdef WIN32
    strcpy(pathBuffer, path);
#endif
}
//------------------------------------------------------------------------------
