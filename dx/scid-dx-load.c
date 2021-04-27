/* Automatically generated - may need to edit! */

#include <dx/dx.h>
#include <dx/modflags.h>

#if defined(intelnt) || defined(WIN32)
#include <windows.h>
#endif

#if defined(__cplusplus)
extern "C" Error DXAddModule (char *, ...);
#else
extern Error DXAddModule (char *, ...);
#endif

#if defined(__cplusplus)
extern "C" Error m_SCID(Object*, Object*);
#endif
#if defined(intelnt) || defined(WIN32)
void FAR WINAPI DXEntry()
#else
  #if defined(__cplusplus)
    extern "C" void DXEntry()
  #else
    void DXEntry()
  #endif
#endif
{
    {
#if defined(__cplusplus)
        extern "C" Error m_SCID(Object *, Object *);
#else
        extern Error m_SCID(Object *, Object *);
#endif
        DXAddModule("SCID", m_SCID, 
            MODULE_PERSISTENT,
            2, "FileName", "PlotGrid",
            1, "PlotData");
    }
    {
#if defined(__cplusplus)
        extern "C" Error m_FLATTEN(Object *, Object *);
#else
        extern Error m_FLATTEN(Object *, Object *);
#endif
        DXAddModule("FLATTEN", m_FLATTEN, 
            MODULE_PERSISTENT,
            3, "Field", "Axis", "Position",
            1, "FlatField");
    }
    {
#if defined(__cplusplus)
        extern "C" Error m_NORM1(Object *, Object *);
#else
        extern Error m_NORM1(Object *, Object *);
#endif
        DXAddModule("NORM1", m_NORM1, 
            MODULE_PERSISTENT,
            1, "Field",
            1, "Norm");
    }
    {
#if defined(__cplusplus)
        extern "C" Error m_UNWRAP(Object *, Object *);
#else
        extern Error m_UNWRAP(Object *, Object *);
#endif
        DXAddModule("UNWRAP", m_UNWRAP, 
            MODULE_PERSISTENT,
            3, "Phase", "Magnitude", "Fancy",
            1, "UnwrapPhase");
    }
    {
#ifndef __cplusplus
        extern Error m_WRAP(Object *, Object *);
#endif
        DXAddModule("WRAP", m_WRAP, 
            MODULE_PERSISTENT,
            1, "Phase",
            1, "WrapPhase");
    }
}
