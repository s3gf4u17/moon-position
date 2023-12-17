#include <vector>
#include <math.h>
#include <iostream>

// 2023,4,17,5,0,0
double centuries(int year, int month, int d, int hour, int minute, int second) {
    double day = (double)d+(double)hour/24.0+(double)minute/24.0/60.0+(double)second/24.0/60.0/60.0;
	if (month<=2) {month+=12;year--;}
	int A = year/100;
	int B = 2-A+A/4;
	if (year<1582||(year==1582&&month<10)||(year==1582&&month==10&day<5.0)) B = 0;
	double jde = (int)(365.25*(year+4716))+(int)(30.6001*(month+1))+day+B-1524.5;
	return (jde-2451545.0)/36525.0;
}

double gmst(double T) { // translate julian centuries to greenich median solar time
    double Tu = T * 36525.0;
    return 280.46061837+360.98564736629*(Tu)+0.000387933*pow(T,2)-pow(T,3)/38710000.0;
}

double range_value(double value) {
	while(value<0.0||value>=360.0) {
		if (value>=360.0) {
			value-=360.0;
		} else {
			value+=360.0;
		}
	}
	return value;
}

double to_rad(double deg) {
    return deg/180.0*M_PI;
}

std::vector<std::vector<int>> table47a_a = {
	{0,0,1,0},
    {2,0,-1,0},
    {2,0,0,0},
    {0,0,2,0},
    {0,1,0,0},
    {0,0,0,2},
    {2,0,-2,0},
    {2,-1,-1,0},
    {2,0,1,0},
    {2,-1,0,0},
    {0,1,-1,0},
    {1,0,0,0},
    {0,1,1,0},
    {2,0,0,-2},
    {0,0,1,2},
    {0,0,1,-2},
    {4,0,-1,0},
    {0,0,3,0},
    {4,0,-2,0},
    {2,1,-1,0},
    {2,1,0,0},
    {1,0,-1,0},
    {1,1,0,0},
    {2,-1,1,0},
    {2,0,2,0},
    {4,0,0,0},
    {2,0,-3,0},
    {0,1,-2,0},
    {2,0,-1,2},
    {2,-1,-2,0},
    {1,0,1,0},
    {2,-2,0,0},
    {0,1,2,0},
    {0,2,0,0},
    {2,-2,-1,0},
    {2,0,1,-2},
    {2,0,0,2},
    {4,-1,-1,0},
    {0,0,2,2},
    {3,0,-1,0},
    {2,1,1,0},
    {4,-1,-2,0},
    {0,2,-1,0},
    {2,2,-1,0},
    {2,1,-2,0},
    {2,-1,0,-2},
    {4,0,1,0},
    {0,0,4,0},
    {4,-1,0,0},
    {1,0,-2,0},
    {2,1,0,-2},
    {0,0,2,-2},
    {1,1,1,0},
    {3,0,-2,0},
    {4,0,-3,0},
    {2,-1,2,0},
    {0,2,1,0},
    {1,1,-1,0},
    {2,0,3,0},
    {2,0,-1,-2}
};

std::vector<std::vector<int>> table47b_a = {
	{0,0,0,1},
    {0,0,1,1},
    {0,0,1,-1},
    {2,0,0,-1},
    {2,0,-1,1},
    {2,0,-1,-1},
    {2,0,0,1},
    {0,0,2,1},
    {2,0,1,-1},
    {0,0,2,-1},
    {2,-1,0,-1},
    {2,0,-2,-1},
    {2,0,1,1},
    {2,1,0,-1},
    {2,-1,-1,1},
    {2,-1,0,1},
    {2,-1,-1,-1},
    {0,1,-1,-1},
    {4,0,-1,-1},
    {0,1,0,1},
    {0,0,0,3},
    {0,1,-1,1},
    {1,0,0,1},
    {0,1,1,1},
    {0,1,1,-1},
    {0,1,0,-1},
    {1,0,0,-1},
    {0,0,3,1},
    {4,0,0,-1},
    {4,0,-1,1},
    {0,0,1,-3},
    {4,0,-2,1},
    {2,0,0,-3},
    {2,0,2,-1},
    {2,-1,1,-1},
    {2,0,-2,1},
    {0,0,3,-1},
    {2,0,2,1},
    {2,0,-3,-1},
    {2,1,-1,1},
    {2,1,0,1},
    {4,0,0,1},
    {2,-1,1,1},
    {2,-2,0,-1},
    {0,0,1,3},
    {2,1,1,-1},
    {1,1,0,-1},
    {1,1,0,1},
    {0,1,-2,-1},
    {2,1,-1,-1},
    {1,0,1,1},
    {2,-1,-2,-1},
    {0,1,2,1},
    {4,0,-2,-1},
    {4,-1,-1,-1},
    {1,0,1,-1},
    {4,0,1,-1},
    {1,0,-1,-1},
    {4,-1,0,-1},
    {2,-2,0,1}
};

std::vector<int> L_coefficients = {
	6288774,
    1274027,
    658314,
    213618,
    -185116,
    -114332,
    58793,
    57066,
    53322,
    45758,
    -40923,
    -34720,
    -30383,
    15327,
    -12528,
    10980,
    10675,
    10034,
    8548,
    -7888,
    -6766,
    -5163,
    4987,
    4036,
    3994,
    3861,
    3665,
    -2689,
    -2602,
    2390,
    -2348,
    2236,
    -2120,
    -2069,
    2048,
    -1773,
    -1595,
    1215,
    -1110,
    -892,
    -810,
    759,
    -713,
    -700,
    691,
    596,
    549,
    537,
    520,
    -487,
    -399,
    -381,
    351,
    -340,
    330,
    327,
    -323,
    299,
    294,
    0
};

std::vector<int> R_coefficients = {
	-20905355,
    -3699111,
    -2955968,
    -569925,
    48888,
    -3149,
    246158,
    -152138,
    -170733,
    -204586,
    -129620,
    108743,
    104755,
    10321,
    0,
    79661,
    -34782,
    -23210,
    -21636,
    24208,
    30824,
    -8379,
    -16675,
    -12831,
    -10445,
    -11650,
    14403,
    -7003,
    0,
    10056,
    6322,
    -9884,
    5751,
    0,
    -4950,
    4130,
    0,
    -3958,
    0,
    3258,
    2616,
    -1897,
    -2117,
    2354,
    0,
    0,
    -1423,
    -1117,
    -1571,
    -1739,
    0,
    -4421,
    0,
    0,
    0,
    0,
    1165,
    0,
    0,
    8752
};

std::vector<int> B_coefficients = {
	5128122,
    280602,
    277693,
    173237,
    55413,
    46271,
    32573,
    17198,
    9266,
    8822,
    8216,
    4324,
    4200,
    -3359,
    2463,
    2211,
    2065,
    -1870,
    1828,
    -1794,
    -1749,
    -1565,
    -1491,
    -1475,
    -1410,
    -1344,
    -1335,
    1107,
    1021,
    833,
    777,
    671,
    607,
    596,
    491,
    -451,
    439,
    422,
    421,
    -366,
    -351,
    331,
    315,
    302,
    -283,
    -229,
    223,
    223,
    -220,
    -220,
    -185,
    181,
    -177,
    176,
    166,
    -164,
    132,
    -119,
    115,
    107
};

bool position(double T, double obsPhi, double obsLam, double obsH) {
	double L_PRIME = range_value(218.3164477+481267.88123421*T-0.0015786*T*T); // mean equinox / mean longitude
	double D = range_value(297.8501921+445267.1114034*T-0.0018819*pow(T,2)+pow(T,3)/545868.0-pow(T,4)/113065000.0); // mean moon elongation
	double M = range_value(357.5291092+35999.0502909*T-0.0001536*pow(T,2)+pow(T,3)/24490000.0); // sun mean anomaly
	double M_PRIME = range_value(134.9633964+477198.8675055*T+0.0087414*pow(T,2)+pow(T,3)/69699.0-pow(T,4)/14712000.0); // moon mean anomaly
	double F = range_value(93.2720950+483202.0175233*T-0.0036539*pow(T,2)-pow(T,3)/3526000.0+pow(T,4)/863310000.0); // moon perigee
	double A1 = range_value(119.75+131.849*T); // action of venus
	double A2 = range_value(53.09+479264.290*T); // action of jupyter
	double A3 = range_value(313.45+481266.484*T); // ???
	double E = 1-0.002516*T-0.0000074*pow(T,2); // earth orbit eccentricity
	double sumL = 3958.0*sin(to_rad(A1))+1962.0*sin(to_rad(L_PRIME-F))+318.0*sin(to_rad(A2)); // sum of periodic terms for the longitude of the Moon
	for (int i=0;i<L_coefficients.size();i++) {
		double multiplier = (table47a_a[i][1]!=0) ? pow(E,abs(table47a_a[i][1])) : 1.0;
		sumL+=L_coefficients[i]*multiplier*sin(to_rad(table47a_a[i][0]*D+table47a_a[i][1]*M+table47a_a[i][2]*M_PRIME+table47a_a[i][3]*F));
	}
	double sumR = 0.0; // sum of periodic terms for the distance of the Moon
	for (int i=0;i<R_coefficients.size();i++) {
		double multiplier = (table47a_a[i][1]!=0) ? pow(E,abs(table47a_a[i][1])) : 1.0;
		sumR+=R_coefficients[i]*multiplier*cos(to_rad(table47a_a[i][0]*D+table47a_a[i][1]*M+table47a_a[i][2]*M_PRIME+table47a_a[i][3]*F));
	}
	double sumB = -2235.0*sin(to_rad(L_PRIME))+382.0*sin(to_rad(A3))+175.0*sin(to_rad(A1-F)); // sum of periodic terms for the latitude of the Moon
	sumB += 175.0*sin(to_rad(A1+F)) + 127.0 * sin(to_rad(L_PRIME-M_PRIME))-115.0 * sin(to_rad(L_PRIME+M_PRIME));
	for (int i=0;i<B_coefficients.size();i++) {
		double multiplier = (table47b_a[i][1]!=0) ? pow(E,abs(table47b_a[i][1])) : 1.0;
		sumB+=B_coefficients[i]*multiplier*sin(to_rad(table47b_a[i][0]*D+table47b_a[i][1]*M+table47b_a[i][2]*M_PRIME+table47b_a[i][3]*F));
	}
	double lambda = L_PRIME+sumL/1000000.0;
	double beta = sumB/1000000.0;
	double r = 385000.56*1000.0+sumR;
	double pi = asin(6378.14/r)/M_PI*180.0;
	double alambda = lambda+0.004610;
	double U = T/100.0;
	double e0 = 23.439291-1.300258*U-0.000431*pow(U,2)+0.555347*pow(U,3)-0.014272*pow(U,4)-0.069353*pow(U,5);
	e0 += -0.010847*pow(U,6)+0.001978*pow(U,7)+0.007742*pow(U,8)+0.001608*pow(U,9)+0.000681*pow(U,10);
	double omega = range_value(125.04452-1934.136261*T+0.0020708*pow(T,2)+pow(T,3)/450000.0);
	double L = range_value(280.4665+36000.7698*T); // mean longitude of the sun
	double ed = 0.002556*cos(to_rad(omega))+0.000158*cos(to_rad(2*L))+0.000028*cos(to_rad(2*L_PRIME))-0.000025*cos(to_rad(2*omega));
	double e = ed+e0;
	double alpha = range_value(atan2(sin(to_rad(lambda))*cos(to_rad(e))-tan(to_rad(beta))*sin(to_rad(e)),cos(to_rad(lambda)))/M_PI*180.0);
	double delta = asin(sin(to_rad(beta))*cos(to_rad(e))+cos(to_rad(beta))*sin(to_rad(e))*sin(to_rad(lambda)))/M_PI*180.0;
    std::cout << alpha << "\t" << delta << std::endl;
    // local horizontal values
    double theta = range_value(100.46061837+36000.770053608*T+0.000387933*pow(T,2)-pow(T,3)/38710000);
    double H = range_value(theta+obsLam-alpha);
    std::cout << theta << std::endl;
    double A = range_value(atan2((sin(to_rad(H))),(cos(to_rad(H))*sin(to_rad(obsPhi))-tan(to_rad(delta))*cos(to_rad(obsPhi))))/M_PI*180.0+180.0);
    double h = asin(sin(to_rad(obsPhi)*sin(to_rad(delta))+cos(to_rad(obsPhi))*cos(to_rad(delta))*cos(to_rad(H))))/M_PI*180.0;
    std::cout << A << "\t" << h << std::endl;
    // bool visible = (h>=0.0&&h<=90.0) ? 1 : 0;
    // std::cout << visible << "\t" << h << std::endl;
}

int main() {
    double time = centuries(2023,4,11,1,0,0); // czas w universal time
    double visible = position(time,52.37613,20.93514,8.0);
    // NEED TO FIX: for h != 0 - az and h errors
    // maybe pass to gmst h m s and work with (int)jd
	return 0;
}