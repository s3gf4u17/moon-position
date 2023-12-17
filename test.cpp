#include <iostream>
#include <math.h>
#include <vector>

double to_rad(double deg) {
    return deg/180.0*M_PI;
}

void calculate_position(int Year, int Month, double Day) {
    // to calculate accurately the position of the Moon for a given instant it is necessary
    // to take into account hundreds of periodic terms int the Moon's longitude, latitude
    // and distance. here we shall limit ourselves to the most important periodic terms
    // the accuracy of the results will be approximately 10" longitude and 4" latitude

    // time is expressed in centuries and thus should be taken with sufficient number
    // of decimals - at least 9 since during 0.000000001 century (3s) the Moon moves
    // over an arc of 1.7 arc seconds
    double time = calculate_time(Year,Month,Day);

    // calculate L', D, M, M' and F angles by means of the following expressions
    // in order to avoid working with large angles, reduce them to less than 2pi

    // Moon's mean longitude (mean equinox of the date), including constant term
    // of the effect of light-time (-0".7)
    double L_PRIME = range_value(218.3164477 + 481267.88123421 * time - 0.0015786 * time * time);
    // mean elongation of the moon
    double D = range_value(297.8501921+445267.1114034*time-0.0018819*pow(time,2)+pow(time,3)/545868.0-pow(time,4)/113065000.0);
    // Sun's mean anomaly
    double M = range_value(357.5291092+35999.0502909*time-0.0001536*pow(time,2)+pow(time,3)/24490000.0);
    // Moon's mean anomaly
    double M_PRIME = range_value(134.9633964+477198.8675055*time+0.0087414*pow(time,2)+pow(time,3)/69699.0-pow(time,4)/14712000.0);
    // Moon's argument of latitude (mean distance of the Moon from its ascending node)
    double F = range_value(93.2720950+483202.0175233*time-0.0036539*pow(time,2)-pow(time,3)/3526000.0+pow(time,4)/863310000.0);

    // three further arguments (in degrees) are needed
    double A1 = range_value(119.75+131.849*time); // action of venus
    double A2 = range_value(53.09+479264.290*time); // action of jupiter
    double A3 = range_value(313.45+481266.484*time);
    // eccentricity of earths orbit around the sun
    double E = 1-0.002516*time-0.0000074*pow(time,2);
    // sum of periodic terms for the longitude of the Moon - unit 0.000001 deg
    double SUM_l = 0.0;
    for (int i = 0 ; i < SigmaLArguments.size();i++) {
        double val = SigmaLCoefficients[i]*sin(to_rad(SigmaLArguments[i][0]*D+SigmaLArguments[i][1]*M+SigmaLArguments[i][2]*M_PRIME+SigmaLArguments[i][3]*F));
        if (SigmaLArguments[i][1]!=0) val*=pow(E,abs(SigmaLArguments[i][1]));
        SUM_l += val;
    }
    SUM_l += 3958.0*sin(to_rad(A1))+1962.0*sin(to_rad(L_PRIME-F))+318*sin(to_rad(A2));
    // sum of periodic terms for the distance of the Moon - unit 0.001 km
    double SUM_r = 0.0;
    for (int i = 0 ; i < SigmaLArguments.size();i++) {
        double val = SigmaLRCoefficients[i]*cos(to_rad(SigmaLArguments[i][0]*D+SigmaLArguments[i][1]*M+SigmaLArguments[i][2]*M_PRIME+SigmaLArguments[i][3]*F));
        if (SigmaLArguments[i][1]!=0) val*=pow(E,abs(SigmaLArguments[i][1]));
        SUM_r += val;
    }
    // sum of periodic terms for the latitude of the Moon - unit 0.000001 deg
    double SUM_b = 0.0;
    for (int i = 0 ; i < table47b_a.size();i++) {
        double val = table47b_c[i]*sin(to_rad(table47b_a[i][0]*D+table47b_a[i][1]*M+table47b_a[i][2]*M_PRIME+table47b_a[i][3]*F));
        if (table47b_a[i][1]!=0) val*=pow(E,abs(table47b_a[i][1]));
        SUM_b += val;
    }
    SUM_b += -2235.0*sin(to_rad(L_PRIME))+382.0*sin(to_rad(A3))+175.0*sin(to_rad(A1-F));
    SUM_b += 175.0*sin(to_rad(A1+F)) + 127.0 * sin(to_rad(L_PRIME-M_PRIME))-115.0 * sin(to_rad(L_PRIME+M_PRIME));

    double lambda = L_PRIME+SUM_l/1000000.0;
    double beta = SUM_b/1000000.0;
    double delta = 385000.56+SUM_r/1000.0;
    double pi = asin(6378.14/delta)/M_PI*180.0;

    double apparent_lambda = lambda+0.004610;
    // transformation from ecliptical into equatorial coordinates
    double U = time/100;
    double epsilon0 = 23.439291-1.300258*U-0.000431*pow(U,2)+0.555347*pow(U,3)-0.014272*pow(U,4)-0.069353*pow(U,5);
    epsilon0 += -0.010847*pow(U,6)+0.001978*pow(U,7)+0.007742*pow(U,8)+0.001608*pow(U,9)+0.000681*pow(U,10);
    // longitude of ascending node of the Moon's mean orbit on the eclpitic
    double omega = range_value(125.04452-1934.136261*time+0.0020708*pow(time,2)+pow(time,3)/450000.0);
    // mean longitude of the Sun
    double L = range_value(280.4665+36000.7698*time);
    double delta_epsilon = 0.002556*cos(to_rad(omega))+0.000158*cos(to_rad(2*L))+0.000028*cos(to_rad(2*L_PRIME))-0.000025*cos(to_rad(2*omega));
    double epsilon = delta_epsilon+epsilon0;
    double alpha = atan2(sin(to_rad(lambda))*cos(to_rad(epsilon))-tan(to_rad(beta))*sin(to_rad(epsilon)),cos(to_rad(lambda)))/M_PI*180.0;
    double declination = asin(sin(to_rad(beta))*cos(to_rad(epsilon))+cos(to_rad(beta))*sin(to_rad(epsilon))*sin(to_rad(lambda)))/M_PI*180.0;

    std::cout << "lamb: " << apparent_lambda << std::endl;
    std::cout << "beta: " << beta << std::endl;
    std::cout << "dist: " << delta << std::endl;
    std::cout << "alph: " << alpha << std::endl;
    std::cout << "decl: " << declination << std::endl;
}

int main() {
    calculate_position(2023,4,15.9);

    return 0;
}