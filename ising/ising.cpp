
#include "pch.h"

#include "Ising.h"


using namespace std;
namespace fs = std::filesystem;

#define USE_PARALLEL_FOR 1


#if USE_PARALLEL_FOR
//#define MyPolicy std::execution::par_unseq
#define MyPolicy std::execution::par
#else 
#define MyPolicy std::execution::seq
#endif


int L;
double Tmin;
double Tmax;
double T_delta;
double B;
int equilibration_steps;
int datanum;
char name[256];
bool use_parallel;
//const double T_delta = 0.1;

double getM(double T);
void inline LogM(double T);

ofstream tmp_file("tmp_file.txt", ios::out);

struct TMpair {
public:
    double T;
    double M;
    TMpair(double t, double m) : T(t), M(m) {}
};

std::ostream& operator<<(std::ostream& s, const TMpair& tm) {
    s << tm.T << "\t" << tm.M << "\n";
    return s;
}



int main(void)
{
    
    cout << "\nSimulation of the 2D square lattice Ising model,\n enter linear dimension L: ";
    cin >> L;

    cout << "\nenter the min temperature: ";
    cin >> Tmin;

	cout << "\nenter the max temperature: ";
	cin >> Tmax;

	cout << "\nenter the temperature step: ";
	cin >> T_delta;

    cout << "\nenter strength of external magnetic field: ";
    cin >> B;

	cout << "\nHow many Monte Carlo steps per spin for equilibration?: ";
	cin >> equilibration_steps;

    cout << "\nHow many Monte Carlo steps per spin do you want to run?: ";
    cin >> datanum;

    cout << "\nEnter the output file name: ";
    cin >> name;

    // remove duplicate file
    std::filesystem::remove(name);
    
    //use vector container to store the temperature range
    std::vector<double> tmp_data = {};
    for (double tmp = Tmin; tmp < Tmax; tmp += T_delta) {
        tmp_data.push_back(tmp);
    }


    std::for_each(MyPolicy, tmp_data.begin(), tmp_data.end(), LogM);


    tmp_file.close();
    std::filesystem::copy("tmp_file.txt", name);
    std::filesystem::remove("tmp_file.txt");

}



void LogM(double T) {
    tmp_file << TMpair(T, getM(T));
}

double getM(double T) {
    Ising experiment(L, T, B);

    int count = 0;
	do
	{
		experiment.onestep();     // these configurations are discarded to achieve equilibrium
		count++;
	} while (count < equilibration_steps);

    double M_avg = 0;
	count = 0;
	do
	{
		experiment.onestep();
        M_avg += experiment.M/(L*L);
	} while (++count < datanum);

    return(M_avg/datanum);
}



