/**
Implementation of version.h

\file
\copyright Copyright 2018. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATELASTOPLASTICQPOT3D_VERSION_HPP
#define GMATELASTOPLASTICQPOT3D_VERSION_HPP

#include "version.h"

namespace GMatElastoPlasticQPot3d {

namespace detail {

    inline std::string unquote(const std::string& arg)
    {
        std::string ret = arg;
        ret.erase(std::remove(ret.begin(), ret.end(), '\"'), ret.end());
        return ret;
    }

}

inline std::string version()
{
    return detail::unquote(std::string(QUOTE(GMATELASTOPLASTICQPOT3D_VERSION)));
}

inline std::vector<std::string> version_dependencies()
{
    std::vector<std::string> ret;

    ret.push_back("gmatelastoplasticqpot3d=" + version());
    ret.push_back("gmattensor=" + GMatTensor::version());
    ret.push_back("qpot=" + QPot::version());
    ret.push_back("xtensor=" +
        detail::unquote(std::string(QUOTE(XTENSOR_VERSION_MAJOR))) + "." +
        detail::unquote(std::string(QUOTE(XTENSOR_VERSION_MINOR))) + "." +
        detail::unquote(std::string(QUOTE(XTENSOR_VERSION_PATCH))));

    return ret;
}

} // namespace GMatElastoPlasticQPot3d

#endif
