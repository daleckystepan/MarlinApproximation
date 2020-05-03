#include <iostream>
#include <iomanip>
#include <math.h>
#include <float.h>

const int16_t jd_lut_count = 15;
const uint16_t jd_lut_tll = 1 << jd_lut_count;
const int16_t jd_lut_tll0 = __builtin_clz( jd_lut_tll ) + 1;

float jd_lut_k[jd_lut_count] = {};
float jd_lut_b[jd_lut_count] = {};

float jd_lut_k_inv[jd_lut_count] = {};
float jd_lut_b_inv[jd_lut_count] = {};

float jd_lut_A[jd_lut_count] = {};
float jd_lut_B[jd_lut_count] = {};
float jd_lut_C[jd_lut_count] = {};

float jd_marlin(float junction_cos_theta)
{
    const float neg = junction_cos_theta < 0 ? -1 : 1,
    t = neg * junction_cos_theta,
    asinx =       0.032843707f	
                            + t * (-1.451838349f	
                            + t * ( 29.66153956f	
                            + t * (-131.1123477f	
                            + t * ( 262.8130562f	
                            + t * (-242.7199627f + t * 84.31466202f) )))),	
                      junction_theta = M_PI_2 - neg * asinx;
    return junction_theta;
}

inline float approx(float junction_cos_theta)
{
    const float neg = junction_cos_theta < 0 ? -1 : 1,
                      t = neg * junction_cos_theta;

    const int16_t idx = (t == 0.0f) ? 0 : __builtin_clz(uint16_t((1.0f - t) * jd_lut_tll)) - jd_lut_tll0;

    float junction_theta = t * jd_lut_k[idx]+ jd_lut_b[idx];

    if (neg < 0) junction_theta = M_PI - junction_theta;
              return junction_theta;
    return junction_theta;
}

inline float acos_inv(float x)
{
    const float t = (x+1)*0.5;

    const int16_t idx = (t == 0.0f) ? 0 : __builtin_clz(uint16_t((1.0f - t) * jd_lut_tll)) - jd_lut_tll0;

    return x * jd_lut_k_inv[idx] + jd_lut_b_inv[idx];
}

inline float acos_inv_q(float x)
{
    const float t = (x+1)*0.5;

    const int16_t idx = (t == 0.0f) ? 0 : __builtin_clz(uint16_t((1.0f - t) * jd_lut_tll)) - jd_lut_tll0;

    //std::cerr << idx << " " << x <<  " " << x * x * jd_lut_C[idx] + x * jd_lut_B[idx] + jd_lut_A[idx] << std::endl;

    const float ret =  x * (x * jd_lut_C[idx] + jd_lut_B[idx]) + jd_lut_A[idx];  // Improves precision

    return ret;
}

void printLUT(const char *name, float *lut)
{
    std::cerr << std::setprecision(FLT_DECIMAL_DIG) << name << ": " << std::endl;
    for(int i = 0; i < jd_lut_count; ++i) {
        std::cerr << lut[i] << "f, ";
    }
    std::cerr << std::endl; 
}

void printLUT(const char *name, double *lut)
{
    std::cerr << std::setprecision(FLT_DECIMAL_DIG) << name << ": " << std::endl;
    for(int i = 0; i < jd_lut_count; ++i) {
        std::cerr << lut[i] << "f, ";
    }
    std::cerr << std::endl; 
}

int main()
{
    std::cerr << "Generating LUTs" << std::endl;

    double c_approx = 1.00751317f;
    double c_inv = 0.976635218f;
    double c_inv_q = 1; // Q

    for(int i = 0; i < jd_lut_count-1; ++i)
    {
        // approx
        double x0 = (pow(2,i) - 1)/pow(2,i);
        double y0 = acos(x0)*(i==0?1:c_approx);
        double x1 = 0.5*x0 + 0.5;
        double y1 = acos(x1)*c_approx;
        jd_lut_k[i] = (y0-y1)/(x0-x1);
        jd_lut_b[i] = (y1*x0 - y0*x1)/(x0-x1);

        double x2 = 0.5*x0 + 0.5;

        //acos inv
        x0 = 2*x0 - 1;
        y0 = 1.0/acos(x0)*(i==0?1:c_inv);
        x1 = 2*x1 - 1;
        y1 = 1.0/acos(x1)*c_inv;
        jd_lut_k_inv[i] = (y0-y1)/(x0-x1);
        jd_lut_b_inv[i] = (y1*x0 - y0*x1)/(x0-x1);
        
        // Q
        y0 = 1.0/acos(x0)*(i==0?1:c_inv_q);

        x2 = 2*x2 - 1;
        double y2 = 1.0/acos(x2)*c_inv_q;

        // Center between x0 and x2
        x1 = cos(1.0/y0/2 + 1.0/y2/2);
        y1 = 1.0/acos(x1)*c_inv_q;

        double d = 1.0/((x0-x1)*(x0-x2)*(x2-x1));
        //std::cerr << d << std::endl;
        jd_lut_A[i] = - d*x0*x0*x1*y2 + d*x0*x0*x2*y1 + d*x0*x1*x1*y2 - d*x0*x2*x2*y1 - d*x1*x1*x2*y0 + d*x1*x2*x2*y0;
        jd_lut_B[i] = - d*x0*x0*y1 + d*x0*x0*y2 + d*x1*x1*y0 - d*x1*x1*y2 - d*x2*x2*y0 + d*x2*x2*y1;
        jd_lut_C[i] = d*x0*y1 - d*x0*y2 - d*x1*y0 + d*x1*y2 + d*x2*y0 - d*x2*y1;
        
        std::cerr << i << std::setprecision(5) << "\t" << x0 << "\t" << y0 << "\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2
                                               << "\t" << jd_lut_A[i] << "\t" << jd_lut_B[i] << "\t" << jd_lut_C[i] << std::endl;
        //std::cerr << std::setprecision(FLT_DECIMAL_DIG) << jd_lut_k[i] << "\t" << jd_lut_b[i] << std::endl;
    }
    // Last values for approx
    jd_lut_k[jd_lut_count-1] = 0;
    jd_lut_b[jd_lut_count-1] = 0;

    // Last values for approx-inv
    double x0 = (pow(2,jd_lut_count-1) - 1)/pow(2,jd_lut_count-1);
    x0 = 2*x0 - 1;
    double y0 = 1.0/acos(x0)*c_inv;
    double x1 = 0.999999f; // Maximum according to Marlin condition in planner.cpp
    double y1 = 1.0/acos(x1);
    jd_lut_k_inv[jd_lut_count-1] = (y0-y1)/(x0-x1);
    jd_lut_b_inv[jd_lut_count-1] = (y1*x0 - y0*x1)/(x0-x1);

    // Last values for approx-inv Q
    x0 = (pow(2,jd_lut_count-1) - 1)/pow(2,jd_lut_count-1);
    x0 = 2*x0 - 1;
    y0 = 1.0/acos(x0)*c_inv_q;

    double x2 = 0.999999f;
    double y2 = 1.0/acos(x2);

    // Center between x0 and x2
    x1 = cos(1.0/y0/2 + 1.0/y2/2);
    y1 = 1.0/acos(x1)*c_inv_q;

    double d = 1.0/((x0-x1)*(x0-x2)*(x2-x1));
    jd_lut_A[jd_lut_count-1] = - d*x0*x0*x1*y2 + d*x0*x0*x2*y1 + d*x0*x1*x1*y2 - d*x0*x2*x2*y1 - d*x1*x1*x2*y0 + d*x1*x2*x2*y0;
    jd_lut_B[jd_lut_count-1] = - d*x0*x0*y1 + d*x0*x0*y2 + d*x1*x1*y0 - d*x1*x1*y2 - d*x2*x2*y0 + d*x2*x2*y1;
    jd_lut_C[jd_lut_count-1] = d*x0*y1 - d*x0*y2 - d*x1*y0 + d*x1*y2 + d*x2*y0 - d*x2*y1;

    std::cerr << "L: " << std::setprecision(5) << "\t" << x0 << "\t" << y0 << "\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << std::endl;

    //printLUT("K", jd_lut_k);
    //printLUT("B", jd_lut_b);
    //printLUT("K-inv", jd_lut_k_inv);
    //printLUT("B-inv", jd_lut_b_inv);
    printLUT("A", jd_lut_A);
    printLUT("B", jd_lut_B);
    printLUT("C", jd_lut_C);

    std::cerr << "Plotting data" << std::endl;
    double min = 1.0f/0.0f;
    double max = -min;
    for(double t = 0.65; t < 0.999999f; t += 0.000025)
    {
        std::cout << t << " " 
                  << 1.0/acos(t) << " " 
                  << 1.0f/jd_marlin(t) << " " << 1.0/acos(t) - 1.0f/jd_marlin(t) << " "
                  << 1.0f/approx(t)    << " " << 1.0/acos(t) - 1.0f/approx(t)    << " "
                  << acos_inv(t)      << " " << 1.0/acos(t) - acos_inv(t)      << " "
                  << acos_inv_q(t)    << " " << 1.0/acos(t) - acos_inv_q(t)    << " "
                  << std::endl;

        double e = 1.0/(acos(t)*acos_inv_q(t));
        if( isfinite(e) && t < 0.99 )
        {
            min = std::min(min, e);
            max = std::max(max, e);
        }
    }

    std::cerr << std::setprecision(FLT_DECIMAL_DIG) << "Correction factor: " << (min+max)/2.0f << std::endl;

    return 0;
}