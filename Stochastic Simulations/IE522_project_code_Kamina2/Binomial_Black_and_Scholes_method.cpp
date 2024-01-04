// Binomial_Black_and_Scholes_method.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <ctime>
using namespace std;

double expiration, rf, vol, s0, k, q, dt, u, pu, discount_factor;

int n_state;

double** mapper;//For memoization

double N(double z) {
    if (z > 6.0) { return 1.0; };// this guards against overflow 
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a = fabs(z);
    double t = 1.0 / (1.0 + a * p);
    double b = c2 * exp((-z) * (z / 2.0));
    double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
    n = 1.0 - b * n;
    if (z < 0.0) n = 1.0 - n;
    return n;
}

double option_price_put_black_scholes(
    double s,      // spot price
    double expiration) {
    double time_sqrt = sqrt(expiration);
    double d1 = (log(s / k) + (rf - q) * expiration) / (vol * time_sqrt) + 0.5 * vol * time_sqrt;
    double d2 = d1 - (vol * time_sqrt);
    return k * exp(-rf * expiration) * N(-d2) - s * exp(-q * expiration) * N(-d1);
}

double BBS_for_american_put_option(
    int state_ind,
    int up_tick_count,
    double s) {
    //cout << "mapper[" << state_ind << "][" << up_tick_count + state_ind << "] " << mapper[state_ind][up_tick_count + state_ind] << endl;
    if (mapper[state_ind][up_tick_count + state_ind] != -1.) {
        return mapper[state_ind][up_tick_count + state_ind];
    }
    //cout << "state_ind= " << state_ind << " up_tick_count= " << up_tick_count << " s= " << s << endl;
    double tmp_out;
    if (state_ind == n_state - 1) {
        tmp_out = option_price_put_black_scholes(s, dt);
        //cout << "option_price_put_black_scholes with s=" << s << " k=" << k << " rf=" << rf << " vol=" << vol << " dt=" << dt << " price=" << tmp_out << endl;
        mapper[state_ind][up_tick_count + state_ind] = tmp_out;
        return tmp_out;
    }

    double expected_val_if_wait = 
        discount_factor * ( pu * BBS_for_american_put_option(state_ind + 1, up_tick_count + 1, s * u) +
                            (1. - pu) * BBS_for_american_put_option(state_ind + 1, up_tick_count - 1, s / u)
            );

    tmp_out= max(
        max(0., k - s),
        expected_val_if_wait
    );
    mapper[state_ind][up_tick_count + state_ind] = tmp_out;
    return tmp_out;
}



void reinitialize_mapper(double** mapper) {//Just to reuse mapper again when pricing other new options
    for (int i = 0; i < n_state; i++)
        for (int j = 0; j < 2 * i + 1; j++)
            mapper[i][j] = -1;
}

int main(int argc, char* argv[])
{
    //The commented parameters are from the reference paper below 
    //(I used these parameters to check if the result in line with that of the paper)
    //Broadie, M. and Detemple, J., 1996. American option valuation: new bounds, approximations, and a comparison of existing methods. 
    //The Review of Financial Studies, 9(4), pp.1211-1250.
    expiration = 1. / 12.;//0.5;
    rf = 0.04;//0.05;
    vol = 0.2;//0.3;
    s0 = 100.;//100.;
    k = 100.;//90.;
    q = 0.02;//0.;
    n_state = atoi(argv[1]);
    dt = expiration / (double)n_state;
    u = exp(vol * sqrt(dt));
    pu = (exp((rf - q) * dt) - 1. / u) / (u - 1. / u);
    discount_factor = exp(-rf * dt);
    cout << "n_state: " << n_state << endl;
    cout << "dt= " << dt  << " u= " << u << " pu= " << pu << " discount_factor= " << discount_factor << endl;

    clock_t time_before, time_after, time_before2, time_after2;
    double time_used1, time_used2, time_used3;

    //initialize the 2-d mapper array with -1 for memoization (the option value will never be negative, 
    //so if negative, we have not reached that state yet.)
    
    time_before2 = clock();

    mapper = new double* [n_state+1];
    for (int i = 0; i < n_state; i++)
        mapper[i] = new double[2 * i + 1];
    reinitialize_mapper(mapper);// initialize every elements with -1
    time_before = clock();
    double put_price_BBS_n = BBS_for_american_put_option(0, 0, s0);
    time_after = clock();
    time_used1 = ((double)time_after - (double)time_before) / CLOCKS_PER_SEC;

    cout << "BBS put_price (n=" << n_state << "): " << put_price_BBS_n << " time used: " << time_used1 << endl;
    reinitialize_mapper(mapper);// initialize every elements with -1
    n_state = n_state / 2;
    dt = expiration / (double)n_state;
    u = exp(vol * sqrt(dt));
    pu = (exp((rf - q) * dt) - 1. / u) / (u - 1. / u);
    discount_factor = exp(-rf * dt);
    cout << "n_state: " << n_state << endl;
    cout << "dt= " << dt << " u= " << u << " pu= " << pu << " discount_factor= " << discount_factor << endl;
    time_before = clock();
    double put_price_BBS_n_over_2 = BBS_for_american_put_option(0, 0, s0);
    time_after = clock();
    time_used2 = ((double)time_after - (double)time_before) / CLOCKS_PER_SEC;

    time_after2 = clock();
    time_used3 = ((double)time_after2 - (double)time_before2) / CLOCKS_PER_SEC;
    cout << "BBS put_price (n=" << n_state << "): " << put_price_BBS_n_over_2 << " time used: " << time_used2 << endl;

    cout << "BBSR put_price =  " << 2 * put_price_BBS_n - put_price_BBS_n_over_2 << " total time used: " << time_used3 << endl;
    
    
}