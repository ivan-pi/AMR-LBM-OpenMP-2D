/*
 * amrviz_gl.h - the OpenGL 3.3 core types, constants and function pointer
 * types that amrviz.c uses
 *
 * Declared here instead of taken from <GL/glcorearb.h>, which macOS does not
 * ship (Apple has <OpenGL/gl3.h> without the PFN... typedefs) and which
 * Windows does not ship either. The functions themselves are loaded at run
 * time with SDL_GL_GetProcAddress, so no GL header or library is needed.
 * The values are those of the Khronos registry (glcorearb.h).
 */
#ifndef AMRVIZ_GL_H
#define AMRVIZ_GL_H

#include <stddef.h>

#if defined(_WIN32) && !defined(__CYGWIN__)
#define AMRVIZ_APIENTRY __stdcall
#else
#define AMRVIZ_APIENTRY
#endif

typedef unsigned int  GLenum;
typedef unsigned char GLboolean;
typedef unsigned int  GLbitfield;
typedef int           GLint;
typedef int           GLsizei;
typedef unsigned int  GLuint;
typedef float         GLfloat;
typedef char          GLchar;
typedef unsigned char GLubyte;
typedef ptrdiff_t     GLsizeiptr;

#define GL_FALSE                       0
#define GL_TRIANGLES                   0x0004
#define GL_TRIANGLE_STRIP              0x0005
#define GL_SRC_ALPHA                   0x0302
#define GL_ONE_MINUS_SRC_ALPHA         0x0303
#define GL_BLEND                       0x0BE2
#define GL_UNPACK_ALIGNMENT            0x0CF5
#define GL_TEXTURE_2D                  0x0DE1
#define GL_UNSIGNED_BYTE               0x1401
#define GL_FLOAT                       0x1406
#define GL_RED                         0x1903
#define GL_RGB                         0x1907
#define GL_RENDERER                    0x1F01
#define GL_VERSION                     0x1F02
#define GL_NEAREST                     0x2600
#define GL_LINEAR                      0x2601
#define GL_TEXTURE_MAG_FILTER          0x2800
#define GL_TEXTURE_MIN_FILTER          0x2801
#define GL_TEXTURE_WRAP_S              0x2802
#define GL_TEXTURE_WRAP_T              0x2803
#define GL_COLOR_BUFFER_BIT            0x00004000
#define GL_RGB8                        0x8051
#define GL_CLAMP_TO_EDGE               0x812F
#define GL_R8                          0x8229
#define GL_R32F                        0x822E
#define GL_TEXTURE0                    0x84C0
#define GL_TEXTURE1                    0x84C1
#define GL_TEXTURE2                    0x84C2
#define GL_ARRAY_BUFFER                0x8892
#define GL_STREAM_DRAW                 0x88E0
#define GL_FRAGMENT_SHADER             0x8B30
#define GL_VERTEX_SHADER               0x8B31
#define GL_COMPILE_STATUS              0x8B81
#define GL_LINK_STATUS                 0x8B82
#define GL_TEXTURE_BUFFER              0x8C2A
#define GL_MAX_TEXTURE_BUFFER_SIZE     0x8C2B

typedef void (AMRVIZ_APIENTRY *PFNGLCLEARPROC)(GLbitfield mask);
typedef void (AMRVIZ_APIENTRY *PFNGLCLEARCOLORPROC)(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha);
typedef void (AMRVIZ_APIENTRY *PFNGLVIEWPORTPROC)(GLint x, GLint y, GLsizei width, GLsizei height);
typedef void (AMRVIZ_APIENTRY *PFNGLENABLEPROC)(GLenum cap);
typedef void (AMRVIZ_APIENTRY *PFNGLDISABLEPROC)(GLenum cap);
typedef void (AMRVIZ_APIENTRY *PFNGLBLENDFUNCPROC)(GLenum sfactor, GLenum dfactor);
typedef void (AMRVIZ_APIENTRY *PFNGLGETINTEGERVPROC)(GLenum pname, GLint *data);
typedef const GLubyte *(AMRVIZ_APIENTRY *PFNGLGETSTRINGPROC)(GLenum name);
typedef void (AMRVIZ_APIENTRY *PFNGLPIXELSTOREIPROC)(GLenum pname, GLint param);
typedef void (AMRVIZ_APIENTRY *PFNGLGENTEXTURESPROC)(GLsizei n, GLuint *textures);
typedef void (AMRVIZ_APIENTRY *PFNGLBINDTEXTUREPROC)(GLenum target, GLuint texture);
typedef void (AMRVIZ_APIENTRY *PFNGLACTIVETEXTUREPROC)(GLenum texture);
typedef void (AMRVIZ_APIENTRY *PFNGLTEXIMAGE2DPROC)(GLenum target, GLint level, GLint internalformat,
                                                    GLsizei width, GLsizei height, GLint border,
                                                    GLenum format, GLenum type, const void *pixels);
typedef void (AMRVIZ_APIENTRY *PFNGLTEXPARAMETERIPROC)(GLenum target, GLenum pname, GLint param);
typedef void (AMRVIZ_APIENTRY *PFNGLTEXBUFFERPROC)(GLenum target, GLenum internalformat, GLuint buffer);
typedef void (AMRVIZ_APIENTRY *PFNGLGENBUFFERSPROC)(GLsizei n, GLuint *buffers);
typedef void (AMRVIZ_APIENTRY *PFNGLBINDBUFFERPROC)(GLenum target, GLuint buffer);
typedef void (AMRVIZ_APIENTRY *PFNGLBUFFERDATAPROC)(GLenum target, GLsizeiptr size, const void *data, GLenum usage);
typedef void (AMRVIZ_APIENTRY *PFNGLGENVERTEXARRAYSPROC)(GLsizei n, GLuint *arrays);
typedef void (AMRVIZ_APIENTRY *PFNGLBINDVERTEXARRAYPROC)(GLuint array);
typedef void (AMRVIZ_APIENTRY *PFNGLVERTEXATTRIBPOINTERPROC)(GLuint index, GLint size, GLenum type,
                                                             GLboolean normalized, GLsizei stride,
                                                             const void *pointer);
typedef void (AMRVIZ_APIENTRY *PFNGLENABLEVERTEXATTRIBARRAYPROC)(GLuint index);
typedef void (AMRVIZ_APIENTRY *PFNGLVERTEXATTRIBDIVISORPROC)(GLuint index, GLuint divisor);
typedef void (AMRVIZ_APIENTRY *PFNGLDRAWARRAYSPROC)(GLenum mode, GLint first, GLsizei count);
typedef void (AMRVIZ_APIENTRY *PFNGLDRAWARRAYSINSTANCEDPROC)(GLenum mode, GLint first, GLsizei count,
                                                             GLsizei instancecount);
typedef GLuint (AMRVIZ_APIENTRY *PFNGLCREATESHADERPROC)(GLenum type);
typedef void (AMRVIZ_APIENTRY *PFNGLSHADERSOURCEPROC)(GLuint shader, GLsizei count,
                                                      const GLchar *const *string, const GLint *length);
typedef void (AMRVIZ_APIENTRY *PFNGLCOMPILESHADERPROC)(GLuint shader);
typedef void (AMRVIZ_APIENTRY *PFNGLGETSHADERIVPROC)(GLuint shader, GLenum pname, GLint *params);
typedef void (AMRVIZ_APIENTRY *PFNGLGETSHADERINFOLOGPROC)(GLuint shader, GLsizei bufSize, GLsizei *length,
                                                          GLchar *infoLog);
typedef void (AMRVIZ_APIENTRY *PFNGLDELETESHADERPROC)(GLuint shader);
typedef GLuint (AMRVIZ_APIENTRY *PFNGLCREATEPROGRAMPROC)(void);
typedef void (AMRVIZ_APIENTRY *PFNGLATTACHSHADERPROC)(GLuint program, GLuint shader);
typedef void (AMRVIZ_APIENTRY *PFNGLLINKPROGRAMPROC)(GLuint program);
typedef void (AMRVIZ_APIENTRY *PFNGLGETPROGRAMIVPROC)(GLuint program, GLenum pname, GLint *params);
typedef void (AMRVIZ_APIENTRY *PFNGLGETPROGRAMINFOLOGPROC)(GLuint program, GLsizei bufSize, GLsizei *length,
                                                           GLchar *infoLog);
typedef void (AMRVIZ_APIENTRY *PFNGLUSEPROGRAMPROC)(GLuint program);
typedef GLint (AMRVIZ_APIENTRY *PFNGLGETUNIFORMLOCATIONPROC)(GLuint program, const GLchar *name);
typedef void (AMRVIZ_APIENTRY *PFNGLUNIFORM1IPROC)(GLint location, GLint v0);
typedef void (AMRVIZ_APIENTRY *PFNGLUNIFORM1FPROC)(GLint location, GLfloat v0);
typedef void (AMRVIZ_APIENTRY *PFNGLUNIFORM2IPROC)(GLint location, GLint v0, GLint v1);
typedef void (AMRVIZ_APIENTRY *PFNGLUNIFORM2FPROC)(GLint location, GLfloat v0, GLfloat v1);

#endif /* AMRVIZ_GL_H */
