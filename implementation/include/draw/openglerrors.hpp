#pragma once

#if NDEBUG
#define GL_CALL(x) x
#else
#define GL_CALL(x) x getOpenglErrors(__FILE__, __LINE__);
#endif

void getOpenglErrors(const char* file, const unsigned int line);
