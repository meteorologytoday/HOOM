functions {
    real interpolate(real x, real x0, real y0, real x1, real y1) {
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }
}

data {
    int<lower=0> raw_N;     // raw data lengths (monthly)
    int  dom[12];           // days of month
    int  steps[12];         // how many sub intervals from mid-ith-month to mid-(i+1)th-month

    real raw_T[raw_N];      // Average temperature of each month
    real raw_F[raw_N];      // Average heat fluxes of each month

    real T_std;

    real c_sw;
    real rho_sw;
}

transformed data {
    int<lower=2*12+1> N = raw_N - 2 * 12;
    int years = N / 12;
    int sum_steps = sum(steps);
    int total_steps = years * sum_steps;

    real T[N];
    real true_future_T[N];

    real F[total_steps];

    real dt[12];
    real sub_dt[12];

    real heat_cap_density = c_sw * rho_sw;

    if(raw_N % 12 != 0) {
       reject("raw_N = ", raw_N, " is not multiple of 12");
    }

    for(m in 1:12) {
        dt[m] = 86400.0 * (dom[1 + (m - 1) % 12] + dom[1 + (m+1 - 1) % 12]) / 2.0;
        sub_dt[m] = dt[m] / steps[m];
    }

    // The following variables starts from index [period + 2] (e.g. February of 2nd year)
    // because we are making predictioins of next month. So the compared answer is shifted
    // by 1.
    true_future_T = raw_T[12+2:12+2+N-1];
    T             = raw_T[12+1:12+1+N-1];

    // interpolate fluxes
    for(y in 1:years) {
        int t = sum_steps * (y - 1) + 1; 
        int t_month = 12 + (y - 1) * 12 + 1;
        for(m in 1:12) {
            for(k in 1:steps[m]) {
                F[t] = interpolate(k-0.5, 0.0, raw_F[t_month], steps[m], raw_F[t_month+1]);
                //print("F[", t, "] = ", F[t]);
                t += 1;

            }
            t_month += 1;
        }
    }
                
    print("Heat capacity density = ", heat_cap_density);
    print("dom          = ", dom);
    print("steps        = ", steps);
    print("sum_steps    = ", sum_steps);
    print("total_steps  = ", total_steps);
    print("dt           = ", dt);
    print("sub_dt       = ", sub_dt);
    print("mean(F)      = ", mean(F));
    print("sd(F)        = ", sd(F));
    print("            T[13:24]  = ", T[13:24]);
    print("true_future_T[13:24]  = ", true_future_T[13:24]);
}

parameters {
    real<lower=1> h;
    real Q[12];
}

model{

    real heat_cap;
    real Q_interp[sum_steps];

    int t;
    int t_month;
    int step_in_year;

    
    heat_cap = h * heat_cap_density;
    step_in_year = 1;
    for(m in 1:12) {
        for(k in 1:steps[m]) {
            if (m < 12) {
                Q_interp[step_in_year] = interpolate(k-0.5, 0.0, Q[m], steps[m], Q[m+1]);
            } else {
                Q_interp[step_in_year] = interpolate(k-0.5, 0.0, Q[12], steps[12], Q[1]);
            }
    
            step_in_year += 1;
        }
    }

    t = 1;
    t_month = 1;
    for(y in 1:years) {
        step_in_year = 1;
        for(m in 1:12) {
            real predict_T = T[t_month];
            for(k in 1:steps[m]) {

                predict_T += ( ( F[t] + Q_interp[step_in_year] ) / heat_cap ) * sub_dt[m];
               
                step_in_year += 1;
                t += 1;
            }

            // Likelihood
            (predict_T - true_future_T[t_month]) ~ normal(0, T_std);

            t_month += 1;
        }
    }

}
