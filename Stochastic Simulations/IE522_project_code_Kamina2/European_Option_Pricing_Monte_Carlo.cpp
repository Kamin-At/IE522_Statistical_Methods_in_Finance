// European_Option_Pricing_Monte_Carlo.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <tuple>
#include <ctime>
#include <fstream>


using namespace std;

double N(const double& z) {
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
};

double option_price_put_black_scholes(
    double s0,// spot price
    double k,// Strike (exercise) price,
    double rf,// interest rate
    double q,// countinuously compounded dividend
    double sigma,// volatility
    double t)// Remaining time in year
{
    double time_sqrt = sqrt(t);
    double d1 = (log(s0 / k) + (rf - q) * t) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
    double d2 = d1 - (sigma * time_sqrt);
    return k * exp(-rf * t) * N(-d2) - s0 * exp(-q * t) * N(-d1);
};

tuple<double, double, double, double> get_monte_carlo_european_put_option_price_with_variance_reduction(
    double s0,
    double k,
    double t,
    double rf,
    double q,
    double sigma,
    int sample_size) {

    default_random_engine generator;
    normal_distribution<double> distribution(0.0, sigma * sqrt(t));

    double drift_term = (rf - q - 0.5 * pow(sigma, 2.)) * t;
    double risk_free_discount_factor = exp(-rf * t);
    double k_to_s0 = k / s0;
    double option_ratio_value = 0.;
    double sum_square_payoff_value = 0;
    double future_payoff_ratio, future_payoff_ratio2, y;
    double rand_z;
    for (int i = 0; i < sample_size; i++) {
        rand_z = distribution(generator);
        future_payoff_ratio = max(k_to_s0 - exp(drift_term + rand_z), 0.);//European put opotion formula
        future_payoff_ratio2 = max(k_to_s0 - exp(drift_term - rand_z), 0.);//European put opotion formula with negative correlation

        y = (future_payoff_ratio + future_payoff_ratio2) / 2;

        option_ratio_value = (i * option_ratio_value + y) / ((double) i + 1.);
        sum_square_payoff_value = (i * sum_square_payoff_value + pow(y, 2.)) / ((double) i + 1.);

    }



    double standard_error = s0 * risk_free_discount_factor * sqrt((sum_square_payoff_value -
        pow(option_ratio_value, 2)) / ((double)sample_size - 1));
    double option_value = s0 * risk_free_discount_factor * option_ratio_value;

    double upper_bound = option_value + 1.96 * standard_error;//1.96 is used, assuming n is large, so Z approaches 1.96
    double lower_bound = option_value - 1.96 * standard_error;//

    return make_tuple(option_value, standard_error, upper_bound, lower_bound);

}

tuple<double, double, double, double> get_monte_carlo_european_put_option_price(
    double s0,                                    
    double k,
    double t,
    double rf,
    double q,
    double sigma,
    int sample_size) {

    default_random_engine generator;
    normal_distribution<double> distribution(0.0, sigma * sqrt(t));

    double drift_term = (rf - q - 0.5 * pow(sigma, 2.)) * t;
    double risk_free_discount_factor = exp(-rf * t);
    double k_to_s0 = k/s0;
    double option_ratio_value = 0.;
    double sum_square_payoff_value = 0;
    double future_payoff_ratio;
    for (int i = 0; i < sample_size; i++) {
        future_payoff_ratio =  max(k_to_s0 - exp(drift_term + distribution(generator)), 0.);//European put opotion formula
        option_ratio_value = (i * option_ratio_value + future_payoff_ratio) / (i + 1);
        sum_square_payoff_value = (i * sum_square_payoff_value + pow(future_payoff_ratio, 2.)) / (i + 1);
    }
    

    double standard_error = s0 * risk_free_discount_factor * sqrt((sum_square_payoff_value  -
                                 pow(option_ratio_value, 2)) / ((double) sample_size - 1));
    double option_value = s0 * risk_free_discount_factor * option_ratio_value;
    /*cout << "option_value: " << option_value << endl;*/
    double upper_bound = option_value + 1.96 * standard_error;//1.96 is used, assuming n is large, so Z approaches 1.96
    double lower_bound = option_value - 1.96 * standard_error;//

    return make_tuple(option_value, standard_error, upper_bound, lower_bound);
}


int main(int argc, char* argv[])
{
    double s0 = 100;//stod(argv[1]);
    double k = 100;//stod(argv[2]);
    double t = 0.5;//stod(argv[3]);
    double rf = 0.04;//stod(argv[4]);
    double q = 0.02;//stod(argv[5]);
    double sigma = 0.2;//stod(argv[6]);
    //int sample_size = stoi(argv[7]);

    ofstream fout("computing_time_comparison.csv");
    fout << "n,without_variance_reduction_time,without_variance_reduction_price,without_variance_reduction_standard_error,with_variance_reduction_time,with_variance_reduction_price,with_variance_reduction_standard_error\n";
    clock_t time_before, time_after;

    double BSM_put_value = option_price_put_black_scholes(s0, k, rf, q, sigma, t);
    double diff_variance_reduction, diff;
    
    int n[] = {50000, 100000, 500000, 1000000, 5000000, 10000000, 50000000, 100000000, 500000000 };

    for (int i = 0; i < size(n); i++) {
        time_before = clock();
        tuple<double, double, double, double> out_variance_reduction = get_monte_carlo_european_put_option_price_with_variance_reduction(
            s0,
            k,
            t,
            rf,
            q,
            sigma,
            n[i]);
        time_after = clock();
        diff_variance_reduction = ((double)time_after - (double)time_before) / CLOCKS_PER_SEC;

        double option_value_variance_reduction = get<0>(out_variance_reduction);
        double standard_error_variance_reduction = get<1>(out_variance_reduction);
        double upper_bound_variance_reduction = get<2>(out_variance_reduction);
        double lower_bound_variance_reduction = get<3>(out_variance_reduction);
        cout << "n: " << n[i] << " with variance reduction technique, computing time: " << diff_variance_reduction << " BSM_put_value: " << BSM_put_value << endl;
        cout << "option_value: " << option_value_variance_reduction << " standard_error: " << standard_error_variance_reduction << endl;
        cout << " 95% confidence interval: [" << lower_bound_variance_reduction << ", " << upper_bound_variance_reduction << "] with Absolute Error:" << BSM_put_value - option_value_variance_reduction << endl;
        cout << "----------------" << endl;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////// WITH VARIANCE REDUCTION //////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        time_before = clock();
        tuple<double, double, double, double> out = get_monte_carlo_european_put_option_price(
            s0,
            k,
            t,
            rf,
            q,
            sigma,
            n[i]);
        time_after = clock();
        diff = ((double)time_after - (double)time_before) / CLOCKS_PER_SEC;

        double option_value = get<0>(out);
        double standard_error = get<1>(out);
        double upper_bound = get<2>(out);
        double lower_bound = get<3>(out);

        cout << "n: " << n[i] << " without variance reduction technique, computing time: " << diff << " BSM_put_value: " << BSM_put_value << endl;
        cout << "option_value: " << option_value << " standard_error: " << standard_error << endl;
        cout << " 95% confidence interval: [" << lower_bound << ", " << upper_bound << "] with Absolute Error:" << BSM_put_value - option_value << endl;
        cout << "-------------------------------------------------------" << endl;
        fout << n[i] << "," << diff << "," << option_value << "," << standard_error << "," << diff_variance_reduction << "," << option_value_variance_reduction << "," << standard_error_variance_reduction << "\n";
    }
    fout.close();

}