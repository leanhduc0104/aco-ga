#include<bits/stdc++.h>
using namespace std;

#define int long long
#define ii pair<int,int>
#define fi first
#define se second
#define pb push_back


int32_t main(){
    ios_base::sync_with_stdio(0);
    cin.tie(0);cout.tie(0);
    srand(time(0));
    int n = 25;
    string s = "input";
    s = s + to_string(n);
    s = s + ".txt";
    char* tmp;
    tmp = &s[0];
    freopen(tmp, "w", stdout);
   	cout << n << '\n';
   	for(int i = 1; i <= n; i++){
   		int x, y;
   		x = rand() % 60;
   		y = rand() % 60;
   		cout << x << ' ' << y << '\n';
   	}


}