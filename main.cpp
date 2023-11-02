#include<bits/stdc++.h>
#define MAX 1000000000
#define SIZE 100
#define KA 10
#define KB 1
using namespace std;
typedef long long int ll;

vector<vector<int> >road_map;
vector<pair<double, double> > refer[10010];
// refer values: stdev, time, refer index: num1, num2
vector<pair<vector<int>, double> > sequences[10010];
// comparing sequences: sequence, time, stdev

int process_num = 0, N_inp;
ll nck[SIZE][SIZE];
double P[SIZE][SIZE*101], Q[SIZE][SIZE*101], sumP[SIZE][SIZE*101];

void create_refer() {
    for(int i=1; i<=N_inp; i++) {
        refer[i].push_back({-1,-1});
        for(int j=1; j<=N_inp; j++) {
            refer[i].push_back({-1, -1});
        }
    }
}

void create_roadmap(int n, vector<int>& sequence) {
    if(n == 1) {
        vector<int> tmp;
        tmp.push_back(N_inp);
        int sz = sequence.size();
        for(int i=1; i<sz; i++) {
            tmp.push_back(sequence[i]);
            refer[sequence[i-1]][sequence[i]] = {0,0};
        }
        road_map.push_back(tmp);
    }

    for(int i=min(n-1, (n+1)/2); i>=1; i--) {     // n-1 -> min(n-1, (n+1) / 2)
        sequence.push_back(i);
        create_roadmap(i, sequence);
        sequence.pop_back();
    }
}

void create_sequences() {
    pair<vector<int>, double> tmp = {{}, MAX};

    for(int i=1; i<=N_inp; i++) {
        for(int j=0; j<=N_inp; j++) {
            sequences[i].push_back(tmp);
        }
    }
}

double calculateStandardDeviation(const vector<int>& sequence) {
    vector<double> reciprocals;
    for (int i = 0; i < sequence.size(); i++) {
        reciprocals.push_back(1.0 / sequence[i]);
    }

    double mean = 0.0;
    for (double reciprocal : reciprocals) {
        mean += reciprocal;
    }
    mean /= reciprocals.size();

    double variance = 0.0;
    for (double reciprocal : reciprocals) {
        variance += pow(reciprocal - mean, 2);
    }
    variance /= reciprocals.size();

    double standardDeviation = sqrt(variance);

    return standardDeviation;
}

ll calculate_nck(ll i,ll j) {
    if(nck[i][j] != 0) return nck[i][j];
    if(i==j || j==0) return 1;
    if(i < j) return 0;
    return calculate_nck(i-1, j-1) + calculate_nck(i-1, j);
}

void calculatePQ() {
    for(ll i=1; i<=N_inp; i++) {
        for(ll j=1; j<=i; j++) {
            nck[i][j] = calculate_nck(i,j);
        }
    }

    for(ll i=1; i<=N_inp; i++) {
        Q[i][i] = 1;
        for(ll j=1; j<i; j++) {
            Q[i][j] = nck[i][j] / pow(3, i-1);
            Q[i][i] -= Q[i][j];
        }
    }

    for(ll i=1; i<=N_inp; i++) {
        P[i][1] = Q[i][1];
        for(ll j=1; j<=N_inp*100; j++) {
            for(ll k=2; k<=i; k++) P[i][j] += Q[i][k] * P[k][j-1];
        }
    }

    for(ll i=1; i<=N_inp; i++) {
        for(ll j=1; j<=N_inp*100; j++) {
            sumP[i][j] = sumP[i][j-1] + P[i][j];
        }
    }
}

double calculate_time(vector<int>& sequence) {      // time calculation needs to be modified
    double t = 0;
    double phi[SIZE*101];
    phi[0] = 0;
    for(ll k=1; k<=N_inp*100; k++) {
        phi[k] = 1;
        for(auto &j: sequence) {
            phi[k] *= sumP[j][k];
        }
        t += k*(phi[k] - phi[k-1]);
    }

    /*for(auto &j: sequence) {
        printf("%d ", j);
    }
    printf("\n");
    printf("%f\n", t);*/

    return t;
}

void addSequence(int a, int b, vector<int>& sequence) {
    double stdev = calculateStandardDeviation(sequence);
    double average_time = calculate_time(sequence); // not actually, just assuming

    double weight = KA * stdev + KB * average_time;     // KA, KB: defined
    pair<vector<int>, double> tmp = sequences[a][b];
    if(tmp.second > weight) sequences[a][b] = {sequence, weight};
}

void findSequences(int n, int num_left, int min_num, int a, int b, vector<int>& sequence) {
    if(n == 0 && num_left == 0) {
        addSequence(a,b, sequence);
        return;
    }

    if(n < 0 || num_left <= 0) {
        return;
    }

    for(int i = min_num; i <= n; i++) {
        sequence.push_back(i);
        findSequences(n - i, num_left-1, i, a, b, sequence);
        sequence.pop_back();
    }
}

pair<double, int> choose_best() {
    double best_weight = MAX;
    int best_num = 0;
    int road_map_sz = road_map.size();
    for(int it = 0; it < road_map_sz; it++) {
        double cur_weight = 0;
        int sz = road_map[it].size();
        for(int i=1; i<sz; i++) {
            int a = road_map[it][i-1];
            int b = road_map[it][i];
            cur_weight += sequences[a][b].second;
        }
        if(cur_weight < best_weight) {
            best_weight = cur_weight;
            best_num = it;
        }
    }
    return {best_weight, best_num};
}

void print_best(double best_weight, int best_num) {
    for(auto &j: road_map[best_num]) if(j!=1) printf("%d -> ", j);
    printf("1\n");
    int sz = road_map[best_num].size();
    double total_time = 0, total_stdev = 0;
    for(int i=1; i<sz; i++) {
        int a = road_map[best_num][i-1];
        int b = road_map[best_num][i];

        printf("%d -> %d:\n", a,b);

        printf("Sequence: ");
        for(auto &j: sequences[a][b].first) {
            printf("%d ", j);
        }
        printf("\n");

        double time = calculate_time(sequences[a][b].first);
        double stdev = calculateStandardDeviation(sequences[a][b].first);
        printf("time: %f   stdev: %f\n", time, stdev);

        total_time += time;
        total_stdev += stdev;
    }

    printf("The total time is %f while the total stdev is %f\n", total_time, total_stdev);
    printf("The total weight is: %f\n", best_weight);


    printf("%d %f %f\n", N_inp, total_time, total_stdev);
}

int main() {
    scanf("%d", &N_inp);
    calculatePQ();

    vector<int> emp = {};
    vector<pair<vector<int>, pair<double, double> > > seq_emp = {};
    emp.push_back(N_inp);
    create_refer();

    create_roadmap(N_inp, emp);

    create_sequences();

    emp.pop_back();

    for(int i=1; i<=N_inp; i++) {
        for(int j=1; j<i; j++) {
            if(refer[i][j].first == 0) {
                findSequences(i, j, 1, i, j, emp);
            }
        }
    }

    pair<double, int> tmp = choose_best();

    double best_weight = tmp.first;
    int best_num = tmp.second;

    print_best(best_weight, best_num);
}
