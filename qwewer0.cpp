#include <iostream>
#include <iomanip>
#include <math.h>
#include <float.h>

const int32_t jd_lut_count = 21;
const uint32_t jd_lut_tll = 1 << jd_lut_count;
const int32_t jd_lut_tll0 = __builtin_clzl( jd_lut_tll ) + 1;

float jd_lut_k[jd_lut_count] = {};
float jd_lut_b[jd_lut_count] = {};

inline float approx(float junction_cos_theta)
{
    const float neg = junction_cos_theta < 0 ? -1 : 1,
                      t = neg * junction_cos_theta;
    const int32_t idx = (t == 0.0f) ? 0 : __builtin_clzl(uint32_t((1.0f - t) * jd_lut_tll)) - jd_lut_tll0;
    float junction_theta = t * jd_lut_k[idx]+ jd_lut_b[idx];
    if (neg < 0) junction_theta = M_PI - junction_theta;
              return junction_theta;
    return junction_theta;
}

void printLUT(const char *name, float *lut)
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
    for(int i = 0; i < jd_lut_count-1; ++i)
    {
        // approx
        double x0 = (pow(2,i) - 1)/pow(2,i);
        double y0 = acos(x0)*(i==0?1:c_approx);
        double x1 = 0.5*x0 + 0.5;
        double y1 = acos(x1)*c_approx;
        jd_lut_k[i] = (y0-y1)/(x0-x1);
        jd_lut_b[i] = (y1*x0 - y0*x1)/(x0-x1);
    }
    // Last values for approx
    jd_lut_k[jd_lut_count-1] = 0;
    jd_lut_b[jd_lut_count-1] = 0;

    printLUT("K", jd_lut_k);
    printLUT("B", jd_lut_b);

    std::cerr << "Plotting data" << std::endl;
    for(double t = 0; t < 0.999999f; t += 0.000025)
    {
        std::cout << t << " " 
                  << acos(t) << " " 
                  << approx(t) << " " << acos(t) - approx(t) << " "
                  << std::endl;
    }

    return 0;
}