#pragma once
#include <vector>
#include "Vec3.hpp"

template<class T>
struct Field {
    std::vector<T> d;
    Field() = default;
    explicit Field(int n, T v = T{}) : d(n, v) {}
    T&       operator[](int i)       { return d[i]; }
    const T& operator[](int i) const { return d[i]; }
    int size() const { return (int)d.size(); }
};

using volScalar = Field<double>;
using volVector = Field<Vec3>;
