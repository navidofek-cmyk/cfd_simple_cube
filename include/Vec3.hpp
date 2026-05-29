#pragma once

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
    Vec3() = default;
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
};

inline Vec3  operator+(Vec3 a, Vec3 b){ return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3  operator-(Vec3 a, Vec3 b){ return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline Vec3  operator*(double s, Vec3 a){ return {s*a.x, s*a.y, s*a.z}; }
inline Vec3& operator+=(Vec3& a, Vec3 b){ a.x+=b.x; a.y+=b.y; a.z+=b.z; return a; }

inline double mag2(double a){ return a*a; }
inline double mag2(Vec3 a){ return a.x*a.x + a.y*a.y + a.z*a.z; }
