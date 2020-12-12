/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

*/

#ifndef GMATELASTOPLASTICQPOT3D_CONFIG_H
#define GMATELASTOPLASTICQPOT3D_CONFIG_H

#ifdef GMATELASTOPLASTICQPOT3D_ENABLE_ASSERT

    #define GMATELASTOPLASTICQPOT3D_ASSERT(expr) \
        GMATELASTOPLASTICQPOT3D_ASSERT_IMPL(expr, __FILE__, __LINE__)

    #define GMATELASTOPLASTICQPOT3D_ASSERT_IMPL(expr, file, line) \
        if (!(expr)) { \
            throw std::runtime_error( \
                std::string(file) + ':' + std::to_string(line) + \
                ": assertion failed (" #expr ") \n\t"); \
        }

#else

    #define GMATELASTOPLASTICQPOT3D_ASSERT(expr)

#endif

#define GMATELASTOPLASTICQPOT3D_VERSION_MAJOR 0
#define GMATELASTOPLASTICQPOT3D_VERSION_MINOR 11
#define GMATELASTOPLASTICQPOT3D_VERSION_PATCH 0

#define GMATELASTOPLASTICQPOT3D_VERSION_AT_LEAST(x,y,z) \
    (GMATELASTOPLASTICQPOT3D_VERSION_MAJOR > x || (GMATELASTOPLASTICQPOT3D_VERSION_MAJOR >= x && \
    (GMATELASTOPLASTICQPOT3D_VERSION_MINOR > y || (GMATELASTOPLASTICQPOT3D_VERSION_MINOR >= y && \
                                                   GMATELASTOPLASTICQPOT3D_VERSION_PATCH >= z))))

#define GMATELASTOPLASTICQPOT3D_VERSION(x,y,z) \
    (GMATELASTOPLASTICQPOT3D_VERSION_MAJOR == x && \
     GMATELASTOPLASTICQPOT3D_VERSION_MINOR == y && \
     GMATELASTOPLASTICQPOT3D_VERSION_PATCH == z)

#endif
