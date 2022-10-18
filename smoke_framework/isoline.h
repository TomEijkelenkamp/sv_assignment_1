#ifndef ISOLINE_H
#define ISOLINE_H

#include <QVector2D>

#include <array>
#include <functional>
#include <optional>
#include <vector>


class Isoline
{
public:
    enum class AmbiguousCaseDecider
    {
        Midpoint,
        Asymptotic
    };

private:
    std::vector<QVector2D> m_vertices;
    std::vector<float> m_values;
    size_t const m_DIM;
    float const m_isolineRho;
    float const m_cellSideLength;
    QVector2D const m_vertex0;
    Isoline::AmbiguousCaseDecider const m_ambiguousCaseDecider;


public:
    enum class InterpolationMethod
    {
        Linear,
        None
    };

    Isoline(std::vector<float> const &values,
            size_t const valuesSideSize,
            float const isolineRho,
            float const cellSideLength,
            InterpolationMethod const interpolationMethod,
            AmbiguousCaseDecider const ACG);

    std::vector<QVector2D> vertices() const;
    int checkC(float point,float c) const;
    //float Interpolate(float x1,float y1,float x2,float y2,QString  coordinate) const;
    float Interpolate(float l1,float l2,float f1,float f2,float c) const;

};

#endif // ISOLINE_H
