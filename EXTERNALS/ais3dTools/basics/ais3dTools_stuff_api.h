/***************************************************************************
 *  Description: import/export macros for creating DLLS with Microsoft
 *	compiler. Any exported function needs to be declared with the
 *  appropriate G2O_XXXX_API macro. Also, there must be separate macros
 *  for each DLL (arrrrrgh!!!)
 *
 *  17 Jan 2012
 *  Email: pupilli@cs.bris.ac.uk
 ****************************************************************************/
#ifndef AIS3DTOOLS_STUFF_API_H
#define AIS3DTOOLS_STUFF_API_H

#ifdef _MSC_VER
// We are using a Microsoft compiler:

#ifdef AIS3DTOOLS_SHARED_LIBS
#ifdef basics_EXPORTS
#define AIS3DTOOLS_STUFF_API __declspec(dllexport)
#else
#define AIS3DTOOLS_STUFF_API __declspec(dllimport)
#endif
#else
#define AIS3DTOOLS_STUFF_API
#endif

#else
// Not Microsoft compiler so set empty definition:
#define AIS3DTOOLS_STUFF_API
#endif

#endif // AIS3DTOOLS_STUFF_API_H
