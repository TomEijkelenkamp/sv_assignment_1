#include "isoline.h"

#include <QDebug>

Isoline::Isoline(std::vector<float> const &values,
                 size_t const DIM,
                 float const isolineRho,
                 float const cellSideLength,
                 InterpolationMethod const interpolationMethod,
                 AmbiguousCaseDecider const ambiguousCaseDecider)
:
    m_values(values),
    m_DIM(DIM),
    m_isolineRho(isolineRho),
    m_cellSideLength(cellSideLength),
    m_ambiguousCaseDecider(ambiguousCaseDecider)

{
    float px0;
    float px1;
    float px2;
    float px3;
    float py0;
    float py1;
    float py2;
    float py3;
    float v0;
    float v1;
    float v2;
    float v3;
    float c;
    int caseNo;
    for (size_t x=cellSideLength; x<(m_DIM-cellSideLength); x++)
    {
        for (size_t y=cellSideLength; y<(m_DIM-cellSideLength); y++)
        {
            c=isolineRho;                                                   // Threshold
            px0 = x*cellSideLength;                                         // X-coordinate of point #0 in the cell
            py0 = y*cellSideLength;                                         // Y-coordinate of point #0 in the cell
            px1 = x*cellSideLength + cellSideLength;                        // X-coordinate of point #1 in the cell
            py1 = y*cellSideLength;                                         // Y-coordinate of point #1 in the cell
            px2 = x*cellSideLength + cellSideLength;                        // X-coordinate of point #2 in the cell
            py2 = y*cellSideLength + cellSideLength;                        // Y-coordinate of point #2 in the cell
            px3 = x*cellSideLength;                                         // X-coordinate of point #3 in the cell
            py3 = y*cellSideLength + cellSideLength;                        // Y-coordinate of point #3 in the cell

            v0=m_values[x+y*m_DIM];                                         // value of point #0 in the matrix
            v1=m_values[(x+1)+y*m_DIM];                                     // value of point #1 in the matrix
            v2=m_values[(x+1)+(y+1)*m_DIM];                                 // value of point #2 in the matrix
            v3=m_values[x+(y+1)*m_DIM];                                     // value of point #3 in the matrix

            if (checkC(v0,c)== checkC(v1,c) && checkC(v0,c) == checkC(v2,c) && checkC(v0,c) == checkC(v3,c)){caseNo=0;}             // they are all the same
            if (checkC(v0,c)!= checkC(v1,c) && checkC(v1,c) == checkC(v2,c) && checkC(v1,c) == checkC(v3,c)){caseNo=1;}             // v0 is different than the rest
            if (checkC(v0,c)!= checkC(v1,c) && checkC(v0,c) == checkC(v2,c) && checkC(v0,c) == checkC(v3,c)){caseNo=2;}             // v1 is different than the rest
            if (checkC(v0,c)== checkC(v1,c) && checkC(v0,c) == checkC(v3,c) && checkC(v0,c) != checkC(v2,c)){caseNo=3;}             // v2 is different than the rest
            if (checkC(v0,c)== checkC(v1,c) && checkC(v0,c) == checkC(v2,c) && checkC(v0,c) != checkC(v3,c)){caseNo=4;}             // v3 is different than the rest
            if (checkC(v0,c)== checkC(v1,c) && checkC(v2,c) == checkC(v3,c) && checkC(v0,c) != checkC(v2,c)){caseNo=5;}             // v0,v1 are different from v2,v3
            if (checkC(v0,c)== checkC(v2,c) && checkC(v1,c) == checkC(v3,c) && checkC(v0,c) != checkC(v1,c)){caseNo=6;}             // v0,v2 are different from v1,v3
            if (checkC(v0,c)== checkC(v3,c) && checkC(v1,c) == checkC(v2,c) && checkC(v0,c) != checkC(v1,c)){caseNo=7;}             // v0,v3 are different from v1,v2

            switch (interpolationMethod) {
               case InterpolationMethod::Linear:

                if (caseNo==1){                                                                                                     // v0 is different than the rest
                    m_vertices.push_back(QVector2D( Interpolate(px0,px3,v0,v3,c) , Interpolate(py0,py3,v0,v3,c) ));                 // first point  between p0 and p3
                    m_vertices.push_back(QVector2D( Interpolate(px0,px1,v0,v1,c) , Interpolate(py0,py1,v0,v1,c) ));                 // second point between p0 and p1
                }
                if (caseNo==2){                                                                                                     // v1 is different than the rest
                    m_vertices.push_back(QVector2D( Interpolate(px0,px1,v0,v1,c) , Interpolate(py0,py1,v0,v1,c) ));                 // second point between p0 and p1
                    m_vertices.push_back(QVector2D( Interpolate(px1,px2,v1,v2,c) , Interpolate(py1,py2,v1,v2,c) ));                 // second point between p1 and p2
                 }
                 if (caseNo==3){                                                                                                    // v2 is different than the rest
                     m_vertices.push_back(QVector2D( Interpolate(px1,px2,v1,v2,c) , Interpolate(py1,py2,v1,v2,c) ));                // second point between p1 and p2
                     m_vertices.push_back(QVector2D( Interpolate(px2,px3,v2,v3,c) , Interpolate(py2,py3,v2,v3,c) ));                // second point between p2 and p3
                 }
                 if (caseNo==4){                                                                                                    // v3 is different than the rest
                     m_vertices.push_back(QVector2D( Interpolate(px0,px3,v0,v3,c) , Interpolate(py0,py3,v0,v3,c) ));                // first point  between p0 and p3
                     m_vertices.push_back(QVector2D( Interpolate(px2,px3,v2,v3,c) , Interpolate(py2,py3,v2,v3,c) ));                // second point between p2 and p3
                 }
                 if (caseNo==5){                                                                                                    // v0,v1 are different from v2,v3
                     m_vertices.push_back(QVector2D( Interpolate(px0,px3,v0,v3,c) , Interpolate(py0,py3,v0,v3,c) ));                // first point  between p0 and p3
                     m_vertices.push_back(QVector2D( Interpolate(px1,px2,v1,v2,c) , Interpolate(py1,py2,v1,v2,c) ));                // second point between p1 and p2
                 }
                 if (caseNo==7){                                                                                                    // v0,v3 are different from v1,v2
                     m_vertices.push_back(QVector2D( Interpolate(px0,px1,v0,v1,c) , Interpolate(py0,py1,v0,v1,c) ));                // second point between p0 and p1
                     m_vertices.push_back(QVector2D( Interpolate(px2,px3,v2,v3,c) , Interpolate(py2,py3,v2,v3,c) ));                // second point between p2 and p3
                 }

                 if (caseNo==6){
                     switch (ambiguousCaseDecider) {
                     case AmbiguousCaseDecider::Midpoint:
                         if(checkC(float(v0+v1+v2+v3)/4,c)==1){
                                m_vertices.push_back(QVector2D( Interpolate(px0,px3,v0,v3,c) , Interpolate(py0,py3,v0,v3,c) ));     // first point  between p0 and p3
                                m_vertices.push_back(QVector2D( Interpolate(px0,px1,v0,v1,c) , Interpolate(py0,py1,v0,v1,c) ));     // second point between p0 and p1
                                m_vertices.push_back(QVector2D( Interpolate(px2,px3,v2,v3,c) , Interpolate(py2,py3,v2,v3,c) ));     // second point between p2 and p3
                                m_vertices.push_back(QVector2D( Interpolate(px1,px2,v1,v2,c) , Interpolate(py1,py2,v1,v2,c) ));     // second point between p1 and p2
                          }
                          else{
                                m_vertices.push_back(QVector2D( Interpolate(px0,px3,v0,v3,c) , Interpolate(py0,py3,v0,v3,c) ));     // first point  between p0 and p3
                                m_vertices.push_back(QVector2D( Interpolate(px2,px3,v2,v3,c) , Interpolate(py2,py3,v2,v3,c) ));     // second point between p2 and p3
                                m_vertices.push_back(QVector2D( Interpolate(px0,px1,v0,v1,c) , Interpolate(py0,py1,v0,v1,c) ));     // second point between p0 and p1
                                m_vertices.push_back(QVector2D( Interpolate(px1,px2,v1,v2,c) , Interpolate(py1,py2,v1,v2,c) ));     // second point between p1 and p2
                        }
                        break;
                     case AmbiguousCaseDecider::Asymptotic:
                         if(checkC(float((v0-v1)/(v0+v2-v3-v1)),c)==1){
                                m_vertices.push_back(QVector2D( Interpolate(px0,px3,v0,v3,c) , Interpolate(py0,py3,v0,v3,c) ));     // first point  between p0 and p3
                                m_vertices.push_back(QVector2D( Interpolate(px0,px1,v0,v1,c) , Interpolate(py0,py1,v0,v1,c) ));     // second point between p0 and p1
                                m_vertices.push_back(QVector2D( Interpolate(px2,px3,v2,v3,c) , Interpolate(py2,py3,v2,v3,c) ));     // second point between p2 and p3
                                m_vertices.push_back(QVector2D( Interpolate(px1,px2,v1,v2,c) , Interpolate(py1,py2,v1,v2,c) ));     // second point between p1 and p2
                         }
                         else{
                                m_vertices.push_back(QVector2D( Interpolate(px0,px3,v0,v3,c) , Interpolate(py0,py3,v0,v3,c) ));     // first point  between p0 and p3
                                m_vertices.push_back(QVector2D( Interpolate(px2,px3,v2,v3,c) , Interpolate(py2,py3,v2,v3,c) ));     // second point between p2 and p3
                                m_vertices.push_back(QVector2D( Interpolate(px0,px1,v0,v1,c) , Interpolate(py0,py1,v0,v1,c) ));     // second point between p0 and p1
                                m_vertices.push_back(QVector2D( Interpolate(px1,px2,v1,v2,c) , Interpolate(py1,py2,v1,v2,c) ));     // second point between p1 and p2
                      }
                   }
                 }
                break;
               case InterpolationMethod::None:
                    if (caseNo==1){                                                                     // v0 is different than the rest
                        m_vertices.push_back(QVector2D( float(px0+px3)/2  , float(py0+py3)/2    ));     // first point  between p0 and p3
                        m_vertices.push_back(QVector2D( float(px0+px1)/2  , float(py0+py1)/2    ));     // second point between p0 and p1
                    }
                   if (caseNo==2){                                                                      // v1 is different than the rest
                       m_vertices.push_back(QVector2D(  float(px0+px1)/2  , float(py0+py1)/2    ));     // first point  between p0 and p1
                       m_vertices.push_back(QVector2D(  float(px1+px2)/2  , float(py1+py2)/2    ));     // second point between p1 and p2
                    }
                    if (caseNo==3){                                                                     // v2 is different than the rest
                        m_vertices.push_back(QVector2D( float(px1+px2)/2  , float(py1+py2)/2   ));      // first point  between p1 and p2
                        m_vertices.push_back(QVector2D( float(px2+px3)/2  , float(py2+py3)/2   ));      // second point between p2 and p3
                    }
                    if (caseNo==4){                                                                     // v3 is different than the rest
                        m_vertices.push_back(QVector2D( float(px0+px3)/2  , float(py0+py3)/2   ));      // first point  between p0 and p3
                        m_vertices.push_back(QVector2D( float(px2+px3)/2  , float(py2+py3)/2   ));      // second point between p2 and p3
                    }
                    if (caseNo==5){                                                                     // v0,v1 are different from v2,v3
                        m_vertices.push_back(QVector2D( float(px0+px3)/2   , float(py0+py3)/2   ));     // first point  between p0 and p3
                        m_vertices.push_back(QVector2D( float(px1+px2)/2   , float(py1+py2)/2   ));     // second point between p1 and p2
                    }
                    if (caseNo==7){                                                                     // v0,v3 are different from v1,v2
                        m_vertices.push_back(QVector2D( float(px0+px1)/2  , float(py0+py1)/2   ));      // first point  between p0 and p1
                        m_vertices.push_back(QVector2D( float(px2+px3)/2  , float(py2+py3)/2   ));      // second point between p2 and p3
                    }

                    if (caseNo==6){
                        switch (ambiguousCaseDecider) {
                        case AmbiguousCaseDecider::Midpoint:
                            if(checkC(float(v0+v1+v2+v3)/4,c)==1){
                                m_vertices.push_back(QVector2D( float(px0+px3)/2  , float(py0+py3)/2   ));  // first point  between p0 and p3
                                m_vertices.push_back(QVector2D( float(px0+px1)/2  , float(py0+py1)/2   ));  // second point between p0 and p1
                                m_vertices.push_back(QVector2D( float(px2+px3)/2  , float(py2+py3)/2   ));  // first point  between p2 and p3
                                m_vertices.push_back(QVector2D( float(px1+px2)/2  , float(py1+py2)/2   ));  // second point between p1 and p2
                            }
                            else{
                                m_vertices.push_back(QVector2D( float(px0+px3)/2  , float(py0+py3)/2   ));  // first point  between p0 and p3
                                m_vertices.push_back(QVector2D( float(px2+px3)/2  , float(py2+py3)/2   ));  // second point between p2 and p3
                                m_vertices.push_back(QVector2D( float(px0+px1)/2  , float(py0+py1)/2   ));  // first point  between p0 and p1
                                m_vertices.push_back(QVector2D( float(px1+px2)/2  , float(py1+py2)/2   ));  // second point between p1 and p2
                            }
                            break;
                         case AmbiguousCaseDecider::Asymptotic:
                             if(checkC(float((v0-v1)/(v0+v2-v3-v1)),c)==1){
                                m_vertices.push_back(QVector2D( float(px0+px3)/2  , float(py0+py3)/2   ));  // first point  between p0 and p3
                                m_vertices.push_back(QVector2D( float(px0+px1)/2  , float(py0+py1)/2   ));  // second point between p0 and p1
                                m_vertices.push_back(QVector2D( float(px2+px3)/2  , float(py2+py3)/2   ));  // first point  between p2 and p3
                                m_vertices.push_back(QVector2D( float(px1+px2)/2  , float(py1+py2)/2   ));  // second point between p1 and p2
                             }
                             else{
                                m_vertices.push_back(QVector2D( float(px0+px3)/2  , float(py0+py3)/2   ));  // first point  between p0 and p3
                                m_vertices.push_back(QVector2D( float(px2+px3)/2  , float(py2+py3)/2   ));  // second point between p2 and p3
                                m_vertices.push_back(QVector2D( float(px0+px1)/2  , float(py0+py1)/2   ));  // first point  between p0 and p1
                                m_vertices.push_back(QVector2D( float(px1+px2)/2  , float(py1+py2)/2   ));  // second point between p1 and p2
                            }
                       }
                    }
                }
            }
        }
}

std::vector<QVector2D> Isoline::vertices() const
{
    return m_vertices;
}
float Isoline::Interpolate(float l1,float l2,float f1,float f2,float c) const//
{
    float p;
    p = ((f2-c)*l1 + (c-f1)*l2 ) / (f2-f1);
    return p;
}
/*
float Isoline::Interpolate(float x1,float y1,float x2,float y2,QString  coordinate) const//
{
    float d;
    float p;
    //Define the Distance bewteen the two points
    d=sqrt(pow(x1-x2,2)+pow(y1-y2,2));
    //Interpolate for X or Y
    if (coordinate=="x")
       p = x1 + (x2-x1)/d;
    else
       p = y1 + (y2-y1)/d;
    return p;
}
*/
int Isoline::checkC(float point,float c) const
{
    if (point>c) return 1;
    else return 0;
}
