// Asian_option_with_variance_reduction.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <tuple>
#include <random>
#include <deque>
#include <string>
#include <fstream>
#include <ctime>
using namespace std;

default_random_engine generator;
normal_distribution<double> distribution(0.0, 1.);



double simulate_price_path(
    double sigma_sqrt_dt,
    int m,
    double drift_term,
    double dt,
    double s0
    ) 
{        
    double avg_price = 0.;
    double prev_s = s0;
    double tmp_price = 0;
    for (int i = 1; i <= m; i++) {
        tmp_price = prev_s * exp(drift_term * dt + sigma_sqrt_dt * distribution(generator));
        avg_price = (avg_price * ((double) i - 1) + tmp_price) / ((double) i);
        prev_s = tmp_price;
        //cout << "inside => avg_val: " << avg_price << endl;
    }
    return avg_price;
}

tuple<double, double, double> simulate_price_path_to_find_b(
    double sigma_sqrt_dt,
    int m,
    double drift_term,
    double dt,
    double s0,
    double k
)
{
    double avg_xc = 0.;
    double avg_x = 0.;
    double avg_c = 0.;
    double avg_c_sq = 0.;
    double avg_x_sq = 0.;
    double prev_s;
    double x;
    double c;
    double c_sq;
    double x_sq;
    double xc;
    double tmp;

    for (int j = 1; j <= 10000; j++) {
        x = 0.;
        c = 0.;
        c_sq = 0.;
        x_sq = 0.;
        xc = 0.;
        prev_s = s0;
        for (int i = 1; i <= m; i++) {
            tmp = prev_s * exp(drift_term * (dt) + sigma_sqrt_dt * distribution(generator));
            c = (c * ((double)i - 1) + log(tmp)) / ((double)i);
            x = (x * ((double)i - 1) + tmp) / ((double)i);
            prev_s = tmp;
        }
       
        c = max(exp(c) - k, 0.);
        x = max(x - k, 0.);
        xc = x * c;
        //cout << "x: " << x << " c: " << c << " xc: "<< xc <<endl;

        avg_x = (avg_x * ((double)j - 1) + x) / ((double)j);
        avg_c = (avg_c * ((double)j - 1) + c) / ((double)j);
        avg_xc = (avg_xc * ((double)j - 1) + xc) / ((double)j);
        avg_c_sq = (avg_c_sq * ((double)j - 1) + pow(c, 2.)) / ((double)j);
        avg_x_sq = (avg_x_sq * ((double)j - 1) + pow(x, 2.)) / ((double)j);
    }
    double var_c = (avg_c_sq - pow(avg_c, 2.));
    double var_x = (avg_x_sq - pow(avg_x, 2.));
    double est_b_star = (avg_xc - avg_x * avg_c) / var_c;
    double rho_xc = (avg_xc - avg_x * avg_c) / (sqrt(var_c) * sqrt(var_x));

    cout << "b_star: " << est_b_star << " avg_c: " << avg_c << " rho_xc: " << rho_xc;
    return make_tuple(est_b_star, avg_c, rho_xc);
}

tuple<double, double> simulate_price_path_with_geometric_mean(
    double sigma_sqrt_dt,
    int m,
    double drift_term,
    double dt,
    double s0
)   
{
    double avg_price = 0.;
    double tmp_ret = 0.;
    double avg_ret = 0.;
    double prev_s = s0;
    for (int i = 1; i <= m; i++) {
        tmp_ret = prev_s * exp( drift_term * (dt) + sigma_sqrt_dt * distribution(generator));
        avg_ret = (avg_ret * ((double)i - 1) + log(tmp_ret)) / ((double)i);
        avg_price = (avg_price * ((double)i - 1) + tmp_ret) / ((double)i);
        prev_s = tmp_ret;
    }
    //cout << "Arith ret: " << avg_price << " Geo ret: " << exp(avg_ret) << endl;
    return make_tuple(avg_price, exp(avg_ret));
}

tuple<double, double, double, double, double> asian_option_pricing_with_control_variable(
    double s0,
    double k,
    double t,
    double rf,
    double q,
    double sigma,
    int sample_size,
    int m
) 
{
    // Control variable C comes from the geometrically averaged price
    // So, instead of finding the mean of the discounted payoffs (x) directly,
    // we find the mean of y = x - b(C - E[C]) = E[x - bC] + b E[C], where b_star = Cov(x,C)/Var(C) <= b_star came from minimization of Var[y]
    // To ease our computation, we compute E[x], b, and E[C] first.
    // The Var[y] = Var[x - bC] = Var[x] + b^2 * Var[C] - 2b * Cov(x,C) = Var[x] - Cov(x,C)^2 / Var[C].
    // So, the S.E.[y] = sqrt(Var[x] * (1 - Rho(x,C)^2) / n)


    double dt = t / (double)m;
    double drift_term = (rf - q - 0.5 * pow(sigma, 2.));
    double sigma_sqrt_dt = sqrt(dt) * sigma;
    double discount_factor = exp(- rf * t);

    double tmp_val, tmp_arith_val, tmp_geo_val;


    double avg_squared_val = 0.;
    double avg_val = 0;
    /////////////////////// Finding b star and E[C]////////////////////////////////
    
    tuple<double, double, double> tmp_b_star = simulate_price_path_to_find_b(sigma_sqrt_dt, m, drift_term, dt, s0, k);
    double b_star = get<0>(tmp_b_star);
    double rho = get<2>(tmp_b_star);
    double avg_c = get<1>(tmp_b_star);


    for (int i = 0; i < sample_size; i++) {
        tuple <double, double> tmp_ret = simulate_price_path_with_geometric_mean(sigma_sqrt_dt, m, drift_term, dt, s0);
        tmp_geo_val = max(get<1>(tmp_ret) - k, 0.);
        tmp_arith_val = max(get<0>(tmp_ret) - k, 0.);
        tmp_val = tmp_arith_val - b_star * (tmp_geo_val - avg_c);

        avg_squared_val = ((double)i * avg_squared_val + pow(tmp_val, 2.)) / ((double)i + 1.);
        avg_val = ((double)i * avg_val + tmp_val) / ((double)i + 1.);
    }

    double var_val = pow(discount_factor, 2.) * (avg_squared_val  - pow(avg_val, 2.));

    double avg_value = avg_val * discount_factor;
    double standard_error = sqrt(var_val / (sample_size - 1));

    return make_tuple(
        avg_value,
        standard_error, 
        avg_value + 1.96 * standard_error,
        avg_value - 1.96 * standard_error,
        rho
    );//Assuming that LLN and CLT hold true, we use z = 1.96 for 95% confidence level.
    
}

tuple<double, double, double, double> plain_asian_option_pricing(
    double s0,
    double k,
    double t,
    double rf,
    double q,
    double sigma,
    int sample_size,
    int m
)
{
    double dt = t / (double)m;
    double drift_term = (rf - q - 0.5 * pow(sigma, 2.));
    double sigma_sqrt_dt = sqrt(dt) * sigma;
    double avg_val = 0.;
    double discount_factor = exp(-rf * t);

    double tmp_val;
    double avg_squared_val = 0.;
    for (int i = 0; i < sample_size; i++) {
        tmp_val = discount_factor * max(simulate_price_path(sigma_sqrt_dt, m, drift_term, dt, s0) - k, 0.);//call formula
        avg_val = ((double)i * avg_val + tmp_val) / ((double)i + 1);
        avg_squared_val = ((double)i * avg_squared_val + pow(tmp_val, 2.)) / ((double)i + 1.);
    }

    double standard_error = sqrt((avg_squared_val - pow(avg_val, 2)) / (sample_size - 1));

    return make_tuple(
        avg_val,
        standard_error,
        avg_val + 1.96 * standard_error,
        avg_val - 1.96 * standard_error
    );//Assuming that LLN and CLT hold true, we use z = 1.96 for 95% confidence level.

}

//void initialize_price_paths(int m, int sample_size) {
//    price_paths = new double* [sample_size];
//    for (int i = 0; i < sample_size; i++)
//        price_paths[i] = new double[m];
//}

int main(int argc, char* argv[])
{
    ofstream fout("computing_time_comparison_asina_option.csv");
    fout << "n,wo_control_time,wo_control_price,wo_control_standard_error,w_control_time,w_control_price,w_control_standard_error,w_control_rho\n";

    double s0 = 100;//stod(argv[1]);
    double k = 100;// stod(argv[2]);
    double t = 1.;// stod(argv[3]);
    double rf = 0.1;// stod(argv[4]);
    double q = 0.;// stod(argv[5]);
    double sigma = 0.2;// stod(argv[6]);
    int sample_size = 10000;//stoi(argv[7]);
    const int m = 50;//stoi(argv[8]);

    int n[] = { 50000, 100000, 500000, 1000000, 5000000, 10000000};
    clock_t time_before, time_after;
    double time_normal, time_control;
    for (int i=0; i < size(n); i++) {
        /*cout << "before initializing" << endl;*/
        //initialize_price_paths(m, n[i]);
        /*cout << "after initializing" << endl;*/
        time_before = clock();
        tuple<double, double, double, double> out = plain_asian_option_pricing(s0, k, t, rf, q, sigma, n[i], m);
        time_after = clock();
        time_normal = ((double)time_after - (double)time_before) / CLOCKS_PER_SEC;

        time_before = clock();
        tuple<double, double, double, double, double> out_geo = asian_option_pricing_with_control_variable(s0, k, t, rf, q, sigma, n[i], m);
        time_after = clock();
        time_control = ((double)time_after - (double)time_before) / CLOCKS_PER_SEC;

        cout << "n: " << n[i] << endl;
        cout << "without control variable" << endl;
        cout << "option_value: " << get<0>(out) << " standard_error: " << get<1>(out) << " 95% confidence interval: [" << get<3>(out) << ", " << get<2>(out) << "] execution time: " << time_normal << endl;
        cout << "----------------------" << endl;
        cout << "with control variable" << endl;
        cout << "option_value: " << get<0>(out_geo) << " standard_error: " << get<1>(out_geo) << " 95% confidence interval: [" << get<3>(out_geo) << ", " << get<2>(out_geo) << "] execution time: " << time_control << endl;
        cout << "------------------------------------------------------------------" << endl;
        fout << n[i] << "," << time_normal << "," << get<0>(out) << "," << get<1>(out) << "," << time_control << "," << get<0>(out_geo) << "," << get<1>(out_geo) << "," << get<4>(out_geo) << "\n";
    }
    
}