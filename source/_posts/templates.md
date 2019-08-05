---
title: "模板库"
date: 2019-7-31 11:04:10
categories: 笔记
toc: true
---

1. 算法 (AG)
	1. 网络流 (NF)
		1. 网络流线性规划
			1. 上下界网络流
			1. 最小割
			1. 费用流
		1. Edmond-Karp (EK)
			1. Capacity Scaling
			1. Dinic&SAP
			1. HLPP
	1. 差分约束 (DC)
	1. 通用线性规划 (LP)
		1. 单纯形算法
	1. 分数规划 (FP)
	1. 贪心 (GE)
	1. 分治 (DAC)
		1. 序列分治
			1. [CDQ&~~整体二分&线段树分治~~](#CDQ)
			1. [ ] 快速排序
				1. 快速选择
					1. Median of Medians
		1. 树分治
			1. [点/~~边~~分治](#点分治)
			1. 重链剖分
				1. 树上启发式合并
				1. 长链剖分
	1. 搜索 (SR)
		1. 搜索优化与剪枝
			1. [A* (AS)](#Astar)
			1. 迭代加深搜索 (ID)
		1. 折半搜索 (MIM)
	1. 随机化及近似 (RAN)
		1. ~~爬山~~
		1. [模拟退火](#模拟退火)
		1. ~~遗传算法~~
		1. ~~机器学习基础~~
	1. 离线逆序
		1. [莫队算法](#莫队)
	1. [散列 (HS)](#离散化)
1. 数据结构 (DS)
	1. 栈 (STK) => STL
	1. 队列 (QUE) => STL
	1. [散列表 (HST)](#Hash_Table)
	1. 堆 (HP)
		1. 二叉堆
		1. 可并堆
	1. 二叉查找树 (BST)
		1. 堆树 (THP)
		1. [Treap](#Treap)
		1. 红黑树 (RBT)
		1. [替罪羊树 (SGT)](#替罪羊树)
	1. [树状数组 (BIT)](#树状数组)
	1. [线段树 (SGM)](#线段树)
	1. [并查集 (UFS)](#并查集)
		1. 带权并查集
		1. 路径压缩
		1. 按秩合并
	1. [Sparse Table (ST)](#ST)
	1. [K维树 (KDT)](#KD-Tree)
	1. 动态树
		1. 点/边分治树 (DCT)
		1. [Link-Cut Tree (LCT)](#LCT)
		1. 欧拉回路树 (ETT)
		1. AAA Tree & Top Tree
	1. 动态图
1. 图论 (GT)
	1. 最小生成树
		1. Prim 算法 (PM)
        1. Kruskal 算法 (KS)
			1. Boruvka 算法
		1. 最小树形图
			1. 朱-刘算法
		1. 斯坦纳树
	1. 最短路径
		1. [dijkstra 算法 (DJ)](#Dijkstra)
		1. [Bellman-Ford 算法 (BFD)](#SPFA)
			1. Johnson 算法
		1. Floyd 算法 (FL)
	1. 欧拉路&哈密顿路
	1. 连通性
		1. 点/边双连通分量 (BCC)
		1. 强连通性 (SCC)
		1. 支配树
	1. 匹配、划分与覆盖
		1. KM 算法 (KM)
			1. 交错树
			1. 带花树算法 (BL)
			1. Tutte 矩阵与一般图匹配
		1. 覆盖集与独立集
		1. 稳定婚姻问题与 GS 算法 (GS)
		1. Hall 定理
		1. DAG 路径覆盖
		1. Dilworth 定理
	1. 2-SAT
	1. 虚树
	1. 仙人掌
		1. 圆方树
	1. 弦图与区间图
	1. 图的树分解
	1. 最小割
		1. 最小割树
		1. Stoer-Wagner 算法
	1. 平面图
		1. 平面图对偶图
	1. 网格图
1. 计算几何 (CG)
	1. 几何向量
	1. 二维凸包
		1. 凸包算法
			1. 卷包裹法
			1. 动态凸包
		1. 三维凸包
		1. 半平面交
	1. 旋转卡壳
	1. 三角剖分
		1. V图
	1. 路径规划
1. 代数 (AB)
	1. 微积分基础
		1. Simpson 积分算法
	1. 线性代数
		1. 矩阵基础
			1. [高斯消元](#高斯消元)
			1. 拟阵
			1. Matrix-Tree 定理
			1. 线性递推
	1. 多项式与幂级数
		1. [DFT/FFT (FT)](#FFT)
			1. [NTT](#NTT)
			1. Bluestein 算法
		1. [多项式基本运算](#多项式)
			1. 多项式除法
			1. 多项式基本初等函数
		1. FWT (FWT)
			1. 子集变换
	1. 抽象代数
		1. 置换群
			1. Schreier-Sims 算法
1. 数论 (NT)
	1. 同余和整除
		1. 欧几里得算法
			1. 扩展欧几里得算法
			1. 类欧几里得算法
		1. 欧拉定理
		1. 二次剩余
		1. 原根及离散对数
			1. [BSGS](#BSGS&exBSGS)
		1. lucas 定理
	1. 质数与简单数论函数
		1. 埃氏筛
		1. 欧拉筛
		1. 莫比乌斯反演
		1. 数论函数快速求和
			1. 杜教筛
			1. 洲阁筛
		1. 素性测试
			1. Miller-Robin (MR)
			1. Pollard's Rho 因子分解
1. 组合计数 (CE)
	1. 计数原理
		1. 容斥原理
	1. 计数数列
		1. 斯特林数
		1. 卡特兰数
		1. 伯努利数
	1. 生成函数 (GF)
	1. 杨氏矩阵
	1. Burnside 引理
		1. Polya 定理
1. 博弈论与信息论 (GI)
	1. 博弈基础
		1. 组合游戏
			1. 博弈树与 DAG 模型
			1. Sprague-Grundy 函数 (SG)
			1. Nim (NIM)
				1. Nim 积
			1. 威佐夫博弈
		1. 不平等博弈
			1. 超现实数
		1. 不完全信息博弈
	1. 通信与数据压缩
		1. 校验码
		1. 哈夫曼编码
		1. 游程编码		
1. 形式语言，自动机与串处理 (FAS)
	1. 串处理 (STR)
		1. 模式匹配
			1. KMP 算法 (KMP)
			1. AC 自动机
			1. Shift-And 算法
		1. 字典树 (TRI)
			1. 后缀树
				1. 后缀数组 (SA)
				1. 后缀自动机 (SAM)
		1. Border
			1. Periodicity 引理
		1. 回文串
			1. manacher 算法
			1. 回文自动机 (PAM)
	1. 形式语言
		1. 正则表达式 (RE)
		1. 有限状态自动机 (DFA)
	1. 并行计算

##### CDQ
```c++
#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace std;
const int MAXN = 800005;
struct Query
{
    int x, y, op, val, id, pos;
    bool operator < (const Query &a) const
    {
        return x == a.x ? op < a.op : x < a.x;
    }
} Ask[MAXN], tmp[MAXN];
int n, m, Ans[MAXN];
#define lowbit(_) ((_)&(-_))
struct BIT
{
    int Sum[2000005];
    void add(int x, int c)
    {
        for (int i = x; i <= n; i += lowbit(i))
            Sum[i] += c;
    }
    int Query(int x)
    {
        int ans = 0;
        for (int i = x; i > 0; i -= lowbit(i))
            ans += Sum[i];
        return ans;
    }
}bit;
void add()
{
    int x1, y1, x2, y2;
    scanf ("%d%d%d%d", &x1, &y1, &x2, &y2);
    ++Ans[0];
    Ask[++m].pos = Ans[0], Ask[m].x = x1 - 1, Ask[m].y = y1 - 1, Ask[m].val = 1, Ask[m].op = 1;
    Ask[++m].pos = Ans[0], Ask[m].x = x2    , Ask[m].y = y2    , Ask[m].val = 1, Ask[m].op = 1;
    Ask[++m].pos = Ans[0], Ask[m].x = x1 - 1, Ask[m].y = y2    , Ask[m].val =-1, Ask[m].op = 1;
    Ask[++m].pos = Ans[0], Ask[m].x = x2    , Ask[m].y = y1 - 1, Ask[m].val =-1, Ask[m].op = 1;
}
void CDQ(int l, int r)
{
    if (l == r)
        return;
    int mid = l + r >> 1, l1 = l, l2 = mid + 1;
    for (int i = l; i <= r; i++)
    {
        if (Ask[i].id <= mid && !Ask[i].op)
            bit.add(Ask[i].y, Ask[i].val);
        if (Ask[i].id > mid && Ask[i].op)
            Ans[Ask[i].pos] += Ask[i].val * bit.Query(Ask[i].y);
    }
    for (int i = l; i <= r; i++)
        if (Ask[i].id <= mid && !Ask[i].op)
            bit.add(Ask[i].y, -Ask[i].val);
    for (int i = l; i <= r; i++)
    {
        if (Ask[i].id <= mid)
            tmp[l1++] = Ask[i];
        else tmp[l2++] = Ask[i];
    }
    for (int i = l; i <= r; i++)
        Ask[i] = tmp[i];
    CDQ(l, mid);
    CDQ(mid + 1, r);
    return;
}
int main()
{
    freopen("mokia.in", "r", stdin);
    freopen("mokia.out", "w", stdout);
    int op;
    scanf ("%d%d", &op, &n);
    while (1)
    {
        scanf ("%d", &op);
        if (op == 1)
        {
            m++;
            scanf ("%d%d%d", &Ask[m].x, &Ask[m].y, &Ask[m].val);
        }
        else if (op == 2)
            add();
        else break;
    }
    for (int i = 1; i <= m; i++)
        Ask[i].id = i;
    sort(Ask + 1, Ask + m + 1);
    CDQ(1 ,m);
    for (int i = 1; i <= Ans[0]; i++)
        printf ("%d\n", Ans[i]);
}
```

##### 点分治
```c++
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

inline int read()
{
	int x=0,f=1;char ch=getchar();
	while (ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
	while (ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
	return x*f;
}

const int MAXN = 20005;

struct edge
{
	int END, next;
	int v;
}v[MAXN << 1]; 
int first[MAXN], p;
void add(int a, int b, int c)
{
	v[p].END = b;
	v[p].next = first[a];
	v[p].v = c;
	first[a] = p++;
}
int size[MAXN], Max[MAXN], ans, sum, t[3], d[MAXN];
bool vis[MAXN];
int root;
void Get_root(int rt, int fa)
{
	size[rt] = 1, Max[rt] = 0;
	for (int i = first[rt]; i != -1; i = v[i].next)
	{
		if (!vis[v[i].END] && v[i].END != fa)
		{
			Get_root(v[i].END, rt);
			size[rt] += size[v[i].END];
			Max[rt] = max(Max[rt], size[v[i].END]);
		}
	}
	Max[rt] = max(Max[rt], sum - size[rt]);
	if (Max[rt] < Max[root]) root = rt;
}

void dfs(int rt, int fa)
{
	t[d[rt]]++;
	for (int i = first[rt]; i != -1; i = v[i].next)
	{
		if (!vis[v[i].END] && v[i].END != fa)
		{
			d[v[i].END] = (d[rt] + v[i].v) % 3;
			dfs(v[i].END, rt);
		}
	}
}

int Calc(int rt, int x)
{
	t[0] = t[1] = t[2] = 0;
	d[rt] = x;
	dfs(rt, 0);
	return t[1] * t[2] * 2 + t[0] * t[0];
}

void Solve(int rt)
{
	ans += Calc(rt, 0);
	vis[rt] = 1;
	for (int i = first[rt]; i != -1; i = v[i].next)
	{
		if (!vis[v[i].END])
		{
			ans -= Calc(v[i].END, v[i].v);
			root = 0, sum = size[v[i].END];
			Get_root(v[i].END, 0);
			Solve(root);
		}
	}
}

int gcd(int a, int b)
{
	return b == 0 ? a : gcd(b, a % b);
}

int main()
{
	memset (first, -1, sizeof (first));
	int n = read();
	int a, b, c;
	for (int i = 1; i < n; i++)
	{
		a = read(), b = read(), c = read() % 3;
		add(a, b, c);
		add(b, a, c);
	}
	Max[0] = sum = n;
	Get_root(1, 0);
	Solve(root);
	int G = gcd(ans, n * n);
	printf ("%d/%d\n", ans / G, n * n / G);
}
```

##### Astar
```c++
#include <cstdio>
#include <cstring>
#include <iostream>
#include <queue>
using namespace std;
struct edge {
  int END, v, next;
} ZHENG[10005], FAN[10005];
struct A {
  int END, f, g;
  bool operator<(const A &a) const {
    if (f == a.f) return a.g < g;
    return a.f < f;
  }
};
int first_zheng[1005], first_fan[1005], p;
int dis[1005];
int n;
int s;
void add(int a, int b, int c);
void spfa(int a);
int Astar(int k);
int main() {
  freopen("cowjog.in", "r", stdin);
  freopen("cowjog.out", "w", stdout);
  // cin>>n>>m;
  memset(first_fan, -1, sizeof(first_fan));
  memset(first_zheng, -1, sizeof(first_zheng));
  p = 1;
  int m, a, b, c, K;
  scanf("%d%d%d", &n, &m, &K);
  for (int i = 1; i <= m; i++) {
    scanf("%d%d%d", &a, &b, &c);
    if (a > b)
      add(a, b, c);
    else
      add(b, a, c);
  }
  spfa(1);
  for (int k = 1; k <= K; k++) {
    s = 0;
    cout << Astar(k) << endl;
  }
  // while(1);
}
void add(int a, int b, int c) {
  ZHENG[p].END = b;
  ZHENG[p].v = c;
  ZHENG[p].next = first_zheng[a];
  FAN[p].END = a;
  FAN[p].v = c;
  FAN[p].next = first_fan[b];
  first_zheng[a] = p;
  first_fan[b] = p++;
}
void spfa(int a) {
  memset(dis, 0x3f, sizeof(dis));
  dis[a] = 0;
  bool flag[1005] = {0};
  queue<int> p;
  p.push(a);
  flag[a] = 1;
  while (!p.empty()) {
    int tmp = p.front();
    flag[tmp] = 0;
    p.pop();
    for (int i = first_fan[tmp]; i != -1; i = FAN[i].next) {
      if (dis[FAN[i].END] > dis[tmp] + FAN[i].v) {
        dis[FAN[i].END] = dis[tmp] + FAN[i].v;
        if (!flag[FAN[i].END]) {
          flag[FAN[i].END] = 1;
          p.push(FAN[i].END);
        }
      }
    }
  }
}
int Astar(int k) {
  priority_queue<A> p;
  A s1, s2;
  s1.END = n;
  s1.g = 0;
  s1.f = s1.g + dis[n];
  p.push(s1);
  while (!p.empty()) {
    s2 = p.top();
    p.pop();
    if (s2.END == 1) {
      s++;
      if (s == k) return s2.g;
    }
    for (int i = first_zheng[s2.END]; i != -1; i = ZHENG[i].next) {
      s1.END = ZHENG[i].END;
      s1.g = s2.g + ZHENG[i].v;
      s1.f = s1.g + dis[ZHENG[i].END];
      p.push(s1);
    }
  }
  return -1;
}
```

##### 模拟退火
```c++
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <climits>
#include <ctime>
using namespace std;
const double pi = acos(-1.);
const double eps = 1e-10;
const int seed = time(0);
inline int read()
{
    int x=0,f=1;char ch=getchar();
    while (ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
    while (ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
    return x*f;
}
double Rand()
{
    return 1.0 * rand() / RAND_MAX;
}
struct Point
{
    double x, y, k;
    Point(double _x = 0, double _y = 0): x(_x), y(_y) {
        k = atan2(y, x);
    }
    Point operator - (const Point &_x)
    {
        return Point(x - _x.x, y - _x.y);
    }
    double operator * (const Point &_x)
    {
        return x * _x.y - y * _x.x;
    }
    bool operator < (const Point &_x) const
    {
        return k > _x.k;
    }
}a[10], tmp[10];
int cnt = 0;
double area(Point *x, int n)
{
    double ans = 0;
    cnt++;
    sort(x + 1, x + n + 1);
    Point t(0, 0);
    for (int i = 1; i <= n - 1; i++)
        ans += (x[i] - t) * (x[i + 1] - t);
    ans += (x[n] - t) * (x[1] - t);
    return fabs(ans);
}
int r[10];
bool cmp(const int &a, const int &b)
{
    return a > b;
}
double dis(Point &x, Point &y)
{
    return (x.x - y.x) * (x.x - y.x) + (x.y - y.y) * (x.y - y.y);
}
bool cmp1(Point x, Point y)
{
    double s = (x - tmp[1]) * (y - tmp[1]);
    if (fabs(s) <= eps) return dis(x, tmp[1]) < dis(y, tmp[1]);
    return s > 0;
}
int st[10], top;
double Calc(Point *x, int n)
{
	memcpy(tmp, x, sizeof (tmp));
    int k = 0;
    for (int i = 1; i <= n; i++)
        if (tmp[i].x < tmp[k].x && tmp[i].y < tmp[k].y || !k)
            k = i;
    swap(tmp[1], tmp[k]);
    sort(tmp + 2, tmp + n + 1, cmp1);
    top = 0;
    st[++top] = 1;
    for (int i = 2; i <= n; i++)
    {
        while (top > 1 && (tmp[st[top]] - tmp[st[top - 1]]) * (tmp[i] - tmp[st[top]]) < 0) top--;
        st[++top] = i;
    }
    for (int i = 1; i <= top; i++)
        tmp[i] = tmp[st[i]];
    return area(tmp, top);
}
double T = 1e9, mint = 1e-9, delta = 1 - 1e-2;
int loop = 50;
double w[10];
int main()
{
    // freopen ("1.out", "w", stdout);
    srand(seed);
    int n = read();
    for (int i = 1; i <= n; i++) r[i] = read();
    for (int i = 1; i <= n; i++) w[i] = Rand() * 2 * pi;
    for (int i = 1; i <= n; i++) a[i].x = cos(w[i]) * r[i], a[i].y = sin(w[i]) * r[i];
    double ans = Calc(a, n);
    while (T > mint)
    {
        for (int i = 1; i <= loop; i++)
        {
            int t = rand() % n + 1;
            double tmp_w = w[t];
            w[t] = w[t] + Rand() * 2 * pi * T;
            // if (w[t] > 2 * pi) w[t] -= 2 * pi;
            a[t].x = cos(w[t]) * r[t], a[t].y = sin(w[t]) * r[t];
            double tmp = Calc(a, n);
            // printf ("%.2f\n", (tmp - ans) / T);
            // printf ("%.2f\n", exp((tmp - ans) / T));
            if (tmp > ans || exp((tmp - ans) / T) > Rand()) 
                ans = tmp;
            else
            {
                w[t] = tmp_w;
                a[t].x = cos(w[t]) * r[t], a[t].y = sin(w[t]) * r[t];
            }
        }
        T *= delta;
    }
    printf ("%.10f\n", ans / 2);
    // while (1);
}
```
##### 莫队
```c++
#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace std;
int a[50005], in[50005];
int block = 500;
struct Query
{
	int l, r, ID;
	bool operator < (const Query &b) const 
	{
		return in[l] == in[b.l] ? r > b.r : l < b.l;
	}
}Q[200005];
int cnt[1000005];
int out[200005];
int main()
{
	int n, m, ans = 0;
	scanf("%d", &n);
	for (int i = 1; i <= n; i++)
	{
		scanf("%d", &a[i]);	
		in[i] = (i - 1) / block + 1;
		cnt[a[i]]++;
		if(cnt[a[i]] == 1)
			ans++;
	}
	scanf("%d", &m);
	for (int i = 1; i <= m; i++)
		scanf("%d%d", &Q[i].l, &Q[i].r), Q[i].ID = i;
	sort(Q + 1, Q + m + 1);
	int l = 1, r = n + 1;
	for (int i = 1; i <= m; i++)
	{
		while(l < Q[i].l) {cnt[a[l]]--; if(cnt[a[l]] == 0) ans--; l++;}
		while(l > Q[i].l) {l--; cnt[a[l]]++; if(cnt[a[l]] == 1) ans++;}
		while(r > Q[i].r) {cnt[a[r]]--; if(cnt[a[r]] == 0) ans--; r--;}
		while(r < Q[i].r) {r++; cnt[a[r]]++; if(cnt[a[r]] == 1) ans++;}
		out[Q[i].ID] = ans;
	}
	for (int i = 1; i <= m; i++)
		printf("%d\n", out[i]);
}
```

##### 离散化

```c++
cnt = 0;
for (int i = 1; i <= n; i++) {
  a[i] = read();
  Hash[++cnt] = a[i];
}
sort(Hash + 1, Hash + cnt + 1);
cnt = unique(Hash + 1, Hash + cnt + 1) - Hash - 1;
for (int i = 1; i <= n; i++) {
  a[i] = lower_bound(Hash + 1, Hash + cnt + 1, a[i]) - Hash;
}
```
##### Hash_Table
```c++
const int MAXN = 200000, BASE = 76543;
struct Hash_Table
{
    struct edge
    {
        int Id, next;
        long long Max;
    }v[MAXN];
    int first[BASE + 2], p;
    Hash_Table()
    {
        memset (first, -1, sizeof (first));
        p = 0;
    }
    void clear()
    {
        memset (first, -1, sizeof (first));
        p = 0;
    }
    long long &operator[](int x)
    {
        int H = x % BASE;
        for (int i = first[H]; i != -1; i = v[i].next)
            if (v[i].Id == x)
                return v[i].Max;
        v[p].Id = x;
        v[p].next = first[H];
        first[H] = p++;
        return v[p - 1].Max = -INF;
    }
}Hash;
```

##### Treap
```c++
struct Treap
{
    struct Node
    {
        int s, key, Max_Num, Max_Size;
        data x;
        Node *ch[2];
        Node(data sa)
        {
            ch[0] = ch[1] = NULL;
            s = 1, x = sa, key = rand();
            Max_Num = Max_Size = 0;
        }
#define size(_) ((_) ? (_)->s : 0)
        void Pushup()
        {
            s = size(ch[0]) + size(ch[1]) + 1;
        }
        void Pushdown()
        {
            if (ch[0])
            {
                ch[0]->Max_Num = max(ch[0]->Max_Num, Max_Num);
                ch[0]->Max_Size = max(ch[0]->Max_Size, Max_Size);
                Max_Morale[ch[0]->x.ID] = max(Max_Morale[ch[0]->x.ID], ch[0]->Max_Num);
                Max_Solidarity[ch[0]->x.ID] = max(Max_Solidarity[ch[0]->x.ID], ch[0]->Max_Size);
            }
            if (ch[1])
            {
                ch[1]->Max_Num = max(ch[1]->Max_Num, Max_Num);
                ch[1]->Max_Size = max(ch[1]->Max_Size, Max_Size);
                Max_Morale[ch[1]->x.ID] = max(Max_Morale[ch[1]->x.ID], ch[1]->Max_Num);
                Max_Solidarity[ch[1]->x.ID] = max(Max_Solidarity[ch[1]->x.ID], ch[1]->Max_Size);
            }
            Max_Num = Max_Size = 0;
        }
    } * root;
    Treap()
    {
        root = NULL;
    }
    Node *Merge(Node *A, Node *B)
    {
        if (!A)
            return B;
        if (!B)
            return A;
        if (A->key < B->key)
        {
            A->Pushdown();
            A->ch[1] = Merge(A->ch[1], B);
            A->Pushup();
            return A;
        }
        else
        {
            B->Pushdown();
            B->ch[0] = Merge(A, B->ch[0]);
            B->Pushup();
            return B;
        }
    }
    typedef pair<Node *, Node *> DNode;
    DNode Split(Node *rt, int k)
    {
        if (!rt)
            return DNode(NULL, NULL);
        DNode o;
        rt->Pushdown();
        if (size(rt->ch[0]) >= k)
        {
            o = Split(rt->ch[0], k);
            rt->ch[0] = o.second;
            rt->Pushup();
            o.second = rt;
        }
        else
        {
            o = Split(rt->ch[1], k - size(rt->ch[0]) - 1);
            rt->ch[1] = o.first;
            rt->Pushup();
            o.first = rt;
        }
        return o;
    }
    Node *kth(int k)
    {
        DNode x = Split(root, k - 1);
        DNode y = Split(x.second, 1);
        Node *ans = y.first;
        root = Merge(Merge(x.first, ans), y.second);
        return ans;
    }
    int Rank(Node *rt, data x)
    {
        if (!rt)
            return 0;
        return x <= rt->x ? Rank(rt->ch[0], x) : Rank(rt->ch[1], x) + size(rt->ch[0]) + 1;
    }
    void Insert(data x)
    {
        int k = Rank(root, x);
        if (root)
        {
            root->Max_Size = max(root->Max_Size, size(root));
            root->Max_Num = max(root->Max_Num, x.x);
            Max_Morale[root->x.ID] = max(root->Max_Num, Max_Morale[root->x.ID]);
            Max_Solidarity[root->x.ID] = max(root->Max_Size, Max_Solidarity[root->x.ID]);
        }
        DNode y = Split(root, k);
        Node *n = new Node(x);
        root = Merge(Merge(y.first, n), y.second);
    }
    void remove(data x)
    {
        int k = Rank(root, x);
        DNode a = Split(root, k);
        DNode b = Split(a.second, 1);
        delete b.first;
        root = Merge(a.first, b.second);
    }
};
```

##### 替罪羊树
```c++
namespace Scapegoat_Tree
{
    const int MAXN = 1e6 + 10;
    const double alpha = 0.756;
    struct Node
    {
        Node *ch[2];
        int key, s, cover;
        bool exist;
        void Pushup()
        {
            s = ch[0]->s + ch[1]->s + exist;
            cover = ch[0]->cover + ch[1]->cover + 1;
        }
        bool isBad()
        {
            return ((ch[0]->cover > cover * alpha + 5) || (ch[1]->cover > cover * alpha + 5));
        }
    };
    struct Stree 
    {
        Node mem_poor[MAXN];
        Node *tail, *root, *null, *ls[MAXN];
        Node *bc[MAXN]; int bc_top, top;
        Node *NewNode(int key)
        {
            Node *o = bc_top ? bc[--bc_top] : tail++;
            o->ch[0] = o->ch[1] = null;
            o->s = o->cover = o->exist = 1;
            o->key = key;
            return o;
        }
        void Travel(Node *rt)
        {
            if (rt == null) return;
            Travel(rt->ch[0]);
            if (rt->exist) ls[++top] = rt;
            else bc[bc_top++] = rt;
            Travel(rt->ch[1]);
        }
        Node *Divide(int l, int r)
        {
            if (l > r) return null;
            int mid = (l + r) >> 1;
            Node *o = ls[mid];
            o->ch[0] = Divide(l, mid - 1);
            o->ch[1] = Divide(mid + 1, r);
            o->Pushup();
            return o;
        }
        void ReBuild(Node *&rt)
        {
            top = 0;
            Travel(rt);
            rt = Divide(1, top);
        }
        Node ** Insert(Node *&rt, int val)
        {
            if (rt == null) 
            {
                rt = NewNode(val);
                return &null;
            }
            else
            {
                Node **res = Insert(rt->ch[val >= rt->key], val);
                rt->Pushup();
                if (rt->isBad()) res = &rt;
                return res;
            }
        }
        void erase(Node *rt, int id)
        {
            rt->s--;
            int offset = rt->ch[0]->s + rt->exist;
            if (rt->exist && id == offset)
            {
                rt->exist = false;
                return;
            }
            else
            {
                if (id <= offset) erase(rt->ch[0], id);
                else erase(rt->ch[1], id - offset);
            }
        }
        Stree()
        {
            tail = mem_poor;
            null = tail++;
            null->ch[0] = null->ch[1] = null;
            null->cover = null->s = null->key = 0;
            root = null; bc_top = 0;
        }
        void Insert(int val)
        {
            Node **rt = Insert(root, val);
            if (*rt != null) ReBuild(*rt);
        }
        int Rank(int val)
        {
            Node *rt = root;
            int ans = 1;
            while (rt != null)
            {
                if (rt->key >= val) rt = rt->ch[0];
                else
                {
                    ans += rt->ch[0]->s + rt->exist;
                    rt = rt->ch[1];
                }
            }
            return ans;
        }
        int Kth(int k)
        {
            Node *rt = root;
            while (rt != null)
            {
                if (rt->ch[0]->s + 1 == k && rt->exist) return rt->key;
                else if (rt->ch[0]->s >= k) rt = rt->ch[0];
                else k -= rt->ch[0]->s + rt->exist, rt = rt->ch[1];
            }
        }
        void erase(int k)
        {
            erase(root, Rank(k));
            if (root->s < alpha * root->cover) ReBuild(root);
        }
    };
}
```
##### FFT
```c++
const int MAXN = 60005 << 2;
const double pi = acos(-1.0);
#define Complex complex<double>
int r[MAXN], n;
Complex a[MAXN], b[MAXN];
void DFT(Complex *a, int f)
{
	for (int i = 0; i < n; i++)
		if (r[i] > i) swap(a[i], a[r[i]]);
	for (int i = 1; i < n; i <<= 1)
	{
		Complex wn(cos(pi / i), f * sin(pi / i));
		for (int j = 0; j < n; j += i << 1)
		{
			Complex w = 1;
			for (int k = 0; k < i; k++, w *= wn)
			{
				Complex x = a[j + k], y = w * a[j + k + i];
				a[j + k] = x + y, a[j + k + i] = x - y;
			}
		}
	}
}
int c[MAXN], m;
void FFT(Complex *a, Complex *b)
{
	int l = 0;
	for (n = 1, m <<= 1; n < m; n <<= 1) ++l;
	for (int i = 0; i < n; i++) r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
	DFT(a, 1), DFT(b, 1);
	for (int i = 0; i < n; i++) a[i] *= b[i];
	DFT(a, -1);
	for (int i = 0; i < n; i++) a[i] /= n;
	for (int i = 0; i < m; i++) c[i] = (int)(a[i].real() + 0.5);
}
```

##### NTT
```c++
const int MAXN = 2e6;
const long long MOD = 998244353;
long long A[MAXN], B[MAXN], N;
long long Inv;
int rev[MAXN];
long long pow_mod(long long a, long long b, long long P) {
  long long ans = 1;
  while (b) {
    if (b & 1) ans = ans * a % P;
    b >>= 1;
    a = a * a % P;
  }
  return ans;
}
void FFt(long long *a, int op) {
  long long w, wn, t;
  for (int i = 0; i < N; i++)
    if (i < rev[i]) swap(a[i], a[rev[i]]);
  for (int k = 2; k <= N; k <<= 1) {
    wn = pow_mod(3, op == 1 ? (MOD - 1) / k : MOD - 1 - (MOD - 1) / k, MOD);
    for (int j = 0; j < N; j += k) {
      w = 1;
      for (int i = 0; i < (k >> 1); i++, w = w * wn % MOD) {
        t = a[i + j + (k >> 1)] * w % MOD;
        a[i + j + (k >> 1)] = (a[i + j] - t + MOD) % MOD;
        a[i + j] = (a[i + j] + t) % MOD;
      }
    }
  }
  if (op == -1)
    for (int i = 0; i < N; i++) a[i] = a[i] * Inv % MOD;
}
void FFt(const int *a, const int *b, int *res, int n) {
  N = 1;
  while (N < n) N <<= 1;
  Inv = pow_mod(N, MOD - 2, MOD);
  for (int i = 0; i < N; i++)
    if (i & 1)
      rev[i] = (rev[i >> 1] >> 1) | (N >> 1);
    else
      rev[i] = (rev[i >> 1] >> 1);
  for (int i = 0; i < N; i++) A[i] = a[i], B[i] = b[i];
  FFt(A, 1), FFt(B, 1);
  for (int i = 0; i < N; i++) A[i] = A[i] * B[i] % MOD;
  FFt(A, -1);
  for (int i = 0; i < N; i++) res[i] = A[i];
}
```

##### 树状数组
```c++
long long Sum[MAXN];
int Num[MAXN];
#define lowbit(_) ((_) & (-_))
void add(int x, int c) {
  for (int i = x; i <= cnt; i += lowbit(i)) Sum[i] += c, Num[i]++;
}
long long Query(int x) {
  long long ans = 0;
  for (int i = x; i > 0; i -= lowbit(i))
    ans += Sum[i];
  return ans;
}
int query(int x) {
  int ans = 0;
  for (int i = x; i > 0; i -= lowbit(i))
    ans += Num[i];
  return ans;
}
```

##### 线段树
```c++
const int N = 100005;
int Sum[N << 2];
#define lch l, m, rt << 1
#define rch m + 1, r, rt << 1 | 1
void Pushup(int rt)
{
    Sum[rt] = Sum[rt << 1] + Sum[rt << 1 | 1]; 
}
void Update(int x, int c, int l, int r, int rt)
{
    if (l == r)
    {
        Sum[rt] = c;
        return;
    }
    int m = l + r >> 1;
    if (x <= m) Update(x, c, lch);
    else Update(x, c, rch);
    Pushup(rt);
}
void buildtree(int l, int r, int rt)
{
    if (l == r)
    {
        Sum[rt] = 1;
        return;
    }
    int m = l + r >> 1;
    buildtree(lch);
    buildtree(rch);
    Pushup(rt);
}
int find(int w, int l, int r, int rt)
{
    if (l == r)
        return l;
    int m = l + r >> 1;
    if (Sum[rt << 1] >= w)
        return find(w, lch);
    else
        return find(w - Sum[rt << 1], rch);
}
int Query(int w, int L, int R, int l, int r, int rt)
{
    if (L == l && R == r)
    {
        if (Sum[rt] < w)
            return -Sum[rt];
        return find(w, l, r, rt);
    }
    int m = l + r >> 1;
    if (R <= m)
        return Query(w, L, R, lch);
    else if (L > m)
        return Query(w, L, R, rch);
    else
    {
        int s1 = Query(w, L, m, lch);
        if (s1 > 0) return s1;
        int s2 = Query(w + s1, m + 1, R, rch);
        if (s2 > 0) return s2;
        return s1 + s2;
    }    
}
int Query_Sum(int L, int R, int l, int r, int rt)
{
    if (L == l && R == r)
        return Sum[rt];
    int m = l + r >> 1;
    if (R <= m)
        return Query_Sum(L, R, lch);
    if (L > m)
        return Query_Sum(L, R, rch);
    return Query_Sum(L, m, lch) + Query_Sum(m + 1, R, rch);
}
```

##### 并查集
```c++
int find(int x) {
  if (fa[x] != x) fa[x] = find(fa[x]);
  return fa[x];
}
```


##### ST
```c++
int Max[(50005 << 1) + 200][16], Min[(50005 << 1) + 200][16];

int QueryMax(int l, int r)
{
    if (l > N || r < 0) return -0x3f3f3f3f;
    if (l <= 0) l = 1;
    if (r > N) r = N;
    int k = 0;
    while ((1 << k) <= (r - l + 1)) k++; k--;
    return max(Max[l][k], Max[r - (1 << k) + 1][k]);
}

int QueryMin(int l, int r)
{
    if (l > N || r < 0) return 0x3f3f3f3f;
    if (l <= 0) l = 1;
    if (r > N) r = N;
    int k = 0;
    while ((1 << k) <= (r - l + 1)) k++; k--;
    return min(Min[l][k], Min[r - (1 << k) + 1][k]);
}
int main () {
    for (int i = 1; i <= 15; i++)
        for (int j = 1; j <= N; j++)
        {
            Max[j][i] = max(Max[j][i - 1], Max[j + (1 << (i - 1))][i - 1]);
            Min[j][i] = min(Min[j][i - 1], Min[j + (1 << (i - 1))][i - 1]);
        }
}
```

##### KD-Tree
```c++
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cmath>
  
using namespace std;
  
inline int read()
{
    int x=0,f=1;char ch=getchar();
    while (ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
    while (ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
    return x*f;
}
  
const int INF = 0x3f3f3f3f;
const double alpha = 0.756;
const int MAXN = 5e5 + 5;
  
int now;
  
struct Point
{
    int d[2];
    int &operator[](const int &x)
    {
        return d[x];
    }
    inline bool operator < (const Point &x) const
    {
        return d[now] == x.d[now] ? d[now ^ 1] < x.d[now ^ 1] : d[now] < x.d[now];
    }
}a[MAXN], cur;
  
#define dis(_, __) (\
    int(abs(_[0] - __[0]) + abs(_[1] - __[1]))\
)
  
#define size(_) ((_) ? (_)->s : 0)
  
struct Node
{
    Node *ch[2];
    Point v;
    int s, d;
    int Max[2], Min[2];
    Node(Point x)
    {
        ch[0] = ch[1] = NULL;
        v = x;
        s = 1, d = now;
        Max[0] = Min[0] = x[0];
        Max[1] = Min[1] = x[1];
    }
    Node(){;}
    inline bool operator < (const Node &x) const
    {
        return v < x.v;
    }
    bool IsBad()
    {
        return ((size(ch[0]) > s * alpha) || (size(ch[1]) > s * alpha));
    }
    void Pushup(Node *x)
    {
        if (!x) return;
        for (int i = 0; i <= 1; i++) Min[i] = min(Min[i], x->Min[i]);
        for (int i = 0; i <= 1; i++) Max[i] = max(Max[i], x->Max[i]);
        s += x->s;
    }
    int min_dis()
    {
        int ans = 0;
        ans += max(Min[0] - cur[0], 0) + max(cur[0] - Max[0], 0);
        ans += max(Min[1] - cur[1], 0) + max(cur[1] - Max[1], 0);
        return ans;
    }
}*root;
  
inline void Build(Node *&rt, int l, int r, int d = 0)
{
    if (l > r) return;
    int mid = l + r >> 1;
    now = d;
    nth_element(a + l, a + mid, a + r + 1);
    rt = new Node(a[mid]);
    Build(rt->ch[0], l, mid - 1, d ^ 1);
    Build(rt->ch[1], mid + 1, r, d ^ 1);
    rt->s = 1;
    rt->Pushup(rt->ch[0]);
    rt->Pushup(rt->ch[1]);
}
  
Node **res;
  
inline void Insert(Node *&rt)
{
    if (rt == NULL)
    {
        rt = new Node(cur);
        res = NULL;
        return;
    }
    now = rt->d;
    if (cur < rt->v) Insert(rt->ch[0]);
    else Insert(rt->ch[1]);
    rt->s = 1;
    rt->Pushup(rt->ch[0]);
    rt->Pushup(rt->ch[1]);
    if (rt->IsBad()) res = &rt;
}
  
Node *st[MAXN << 1];
int top = 0;
  
inline void Travel(Node *&rt)
{
    if (!rt) return;
    Travel(rt->ch[0]);
    st[++top] = rt;
    Travel(rt->ch[1]);
}
  
inline int cmp(const Node *x, const Node *y)
{
    return x->v < y->v;
}
  
inline Node *Divide(int l, int r, int d)
{
    if (l > r) return NULL;
    int mid = l + r >> 1;
    now = d;
    nth_element(st + l, st + mid, st + r + 1, cmp);
    Node *rt = st[mid];
    rt->ch[0] = Divide(l, mid - 1, d ^ 1);
    rt->ch[1] = Divide(mid + 1, r, d ^ 1);
    rt->s = 1;
    rt->Max[0] = rt->Min[0] = rt->v[0];
    rt->Max[1] = rt->Min[1] = rt->v[1];
    rt->Pushup(rt->ch[0]);
    rt->Pushup(rt->ch[1]);
}
  
inline void ReBuild(Node *&rt)
{
    top = 0;
    int d = rt->d;
    Travel(rt);
    rt = Divide(1, top, d);
}
  
inline void Insert(Point x)
{
    cur = x;
    Insert(root);
    if (res != NULL)
        ReBuild(*res);
}
  
int Min_ans;
  
inline void Query(Node *rt)
{
    if (!rt) return;
    // if (rt->min_dis() > Min_ans) return;
    Min_ans = min(Min_ans, dis(rt->v, cur));
    int dis_l = rt->ch[0] ? rt->ch[0]->min_dis() : INF;
    int dis_r = rt->ch[1] ? rt->ch[1]->min_dis() : INF;
    if (dis_l < dis_r)
    {
        Query(rt->ch[0]);
        if (dis_r < Min_ans) Query(rt->ch[1]);
    }
    else
    {
        Query(rt->ch[1]);
        if (dis_l < Min_ans) Query(rt->ch[0]);
    }
}
  
inline int Query(Point x)
{
    cur = x;
    Min_ans = INF;
    Query(root);
    return Min_ans;
}
  
int main()
{
    // freopen("1.in", "r", stdin);
    // freopen("2.out", "w", stdout);
    int n, m;
    n = read(), m = read();
    for (int i = 1; i <= n; i++)
        a[i][0] = read(), a[i][1] = read();
    Build(root, 1, n);
    Point x;
    while (m--)
    {
        int t = read();
        if (t == 1)
        {
            x[0] = read(), x[1] = read();
            Insert(x);
        }
        else
        {
            x[0] = read(), x[1] = read();
            printf ("%d\n", Query(x));
        }
    }
}
```

##### LCT
```c++
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

inline int read()
{
	int x=0,f=1;char ch=getchar();
	while (ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
	while (ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
	return x*f;
}

const int MAXN = 10010;

struct Node
{
	Node *f, *ch[2];
	bool IsRoot, Mark;
	Node()
	{
		f = ch[0] = ch[1] = NULL;
		IsRoot = 1, Mark = 0;
	}
}*null, *Q[MAXN];
Node *NewNode()
{
	Node *o = new Node();
	o->f = o->ch[0] = o->ch[1] = null;
	return o;
}
bool son(Node *rt)
{
	return rt->f->ch[1] == rt;
}
inline void Pushflip(Node *rt)
{
	if (rt == null) return;
	rt->Mark ^= 1, swap(rt->ch[0], rt->ch[1]);
}
inline void Pushdown(Node *rt)
{
	if (rt->Mark)
	{
		Pushflip(rt->ch[0]);
		Pushflip(rt->ch[1]);
		rt->Mark = 0;
	}
}
inline void rotate(Node *rt)
{
	if (rt->IsRoot) return;
	Node *fa = rt->f, *Grand = fa->f;
	Pushdown(fa), Pushdown(rt);
	int k = son(rt), kk = son(fa);
	fa->ch[k] = rt->ch[k ^ 1];
	if (rt->ch[k ^ 1] != null) rt->ch[k ^ 1]->f = fa;
	rt->ch[k ^ 1] = fa, fa->f = rt, rt->f = Grand;
	if (!fa->IsRoot) Grand->ch[kk] = rt;
	else fa->IsRoot = 0, rt->IsRoot = 1; 
}
inline void Clear(Node *rt)
{
	if (!rt->IsRoot) Clear(rt->f);
	Pushdown(rt);
}

inline void Splay(Node *rt)
{
	for (Clear(rt); !rt->IsRoot; rotate(rt))
		if (!rt->f->IsRoot)
			rotate(son(rt) == son(rt->f) ? rt->f : rt);
}

inline void Access(Node *rt)
{
	for (Node *pre = null; rt != null; pre = rt, rt = rt->f)
	{
		Splay(rt);
		rt->ch[1]->IsRoot = 1;
		pre->IsRoot = 0;
		rt->ch[1] = pre;
	}
}

inline void MakeRoot(Node *rt)
{
	Access(rt);
	Splay(rt);
	Pushflip(rt);
}

inline void link(Node *a, Node *b)
{
	MakeRoot(a);
	a->f = b;
}

inline void cut(Node *a, Node *b)
{
	MakeRoot(a);
	Access(b);
	Splay(b);
	Pushdown(b);
	b->ch[0]->IsRoot = 1;
	b->ch[0]->f = null;
	b->ch[0] = null;
}

inline Node* find(Node *rt)
{
	if (rt->ch[0] != null) return find(rt->ch[0]);
	return rt;
}

inline bool Judge(Node *a, Node *b)
{
	while (a->f != null) a = a->f;
	while (b->f != null) b = b->f;
	return a == b;
}

int main()
{
	int n = read(), m = read();
	null = new Node();
	null->ch[0] = null->ch[1] = null->f = null, null->IsRoot = null->Mark = 0;
	for (int i = 1; i <= n; i++)
		Q[i] = NewNode();
	int a, b;
	char c[10];
	while (m--)
	{
		scanf ("%s", c);
		if (c[0] == 'C')
		{
			a = read(), b = read();
			link(Q[a], Q[b]);
		}
		else if (c[0] == 'Q')
		{
			a = read(), b = read();
			puts(Judge(Q[a], Q[b]) ? "Yes" : "No");
		}
		else
		{
			a = read(), b = read();
			cut(Q[a], Q[b]);
		}
	}
}

```

##### Dijkstra
```c++
#include <bits/stdc++.h>
const int N = 1e6 + 1;
using namespace std;
struct edge {
  int END, next, v;
} v[N * 2];
int first[N], p;
void add(int a, int b, int c) {
  v[p].END = b;
  v[p].v = c;
  v[p].next = first[a];
  first[a] = p++;
}
int n, m, S;
typedef pair<int, int> T;
int dis[N];
int main() {
  scanf("%d%d%d", &n, &m, &S);
	memset (first, -1, sizeof (first));
  for (int i = 1; i <= m; i++) {
    int x, y, z;
    scanf("%d%d%d", &x, &y, &z);
    add(x, y, z);
  }
  memset(dis, 0x3f, sizeof(dis));
  dis[S] = 0;
  T X;
  X.first = 0, X.second = S;
  priority_queue<T, vector<T>, greater<T> > Q;
  Q.push(X);
  while (!Q.empty()) {
    X = Q.top();
    Q.pop();
    int x = X.second;
    if (dis[x] < X.first) continue;
    for (int i = first[x]; i; i = v[i].next) {
      int y = v[i].END;
      if (dis[y] > dis[x] + v[i].v) {
        dis[y] = dis[x] + v[i].v;
        X.first = dis[y];
        X.second = y;
        Q.push(X);
      }
    }
  }
  for (int i = 1; i <= n; i++) printf("%d ", dis[i]);
  return 0;
}

```

##### SPFA
```c++
int dis[N];
bool flag[N];
queue<int> q;
bool Spfa(int x)
{
	memset(dis, 0x3f, sizeof (dis));
	memset(flag, 0, sizeof (flag));
	dis[x] = 0;
	flag[x] = 1;
	q.push(x);
	while (!q.empty())
	{
		int k = q.front();
		q.pop();
		flag[k] = 0;
		for (int i = first[k]; i != -1; i = v[i].next)
		{
			if (dis[v[i].END] > dis[k] + v[i].v)
			{
				dis[v[i].END] = dis[k] + v[i].v;
				if (!flag[v[i].END])
				{
					flag[v[i].END] = 1;
					q.push(v[i].END);
				}
			}
		}
	}
	return 0;
}
```

##### 多项式
```c++
#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace std;
const int MAXN = 1e6 * 4;
const int MOD  = 998244353;
const int L = (1 << 18) + 1;
const int LM = (1 << 16) + 1;
inline int read()
{
    int x=0,f=1;char ch=getchar();
    while (ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
    while (ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
    return x*f;
}
int N, Inv, rev[MAXN];
int Sum = 0, n, m;
struct data
{
    int e, f, g, id;
}s[60005];
long long pow_mod(long long a, int b)
{
    long long ans = 1;
    while (b)
    {
        if (b & 1)
            ans = ans * a % MOD;
        b >>= 1;
        a = a * a % MOD;
    }
    return ans;
}
int ans[MAXN], cnt;
void insert(int a, int b, int c, int d, int id)
{
    if (b == 0)
    {
        ans[id] = ((1ll * a * Sum) + (1ll * c * n)) % MOD * pow_mod(d, MOD - 2) % MOD;
        return;
    }
    s[++cnt].f = (1ll * b * c - 1ll * a * d) % MOD;
    b = pow_mod(b, MOD - 2);
    s[cnt].e = 1ll * a * b % MOD, s[cnt].g = 1ll * d * b % MOD;
    b = 1ll * b * b % MOD; s[cnt].f = 1ll * s[cnt].f * b % MOD;
    if (s[cnt].f < 0) s[cnt].f += MOD;
    s[cnt].id = id;
}
void Init(int x)
{
    N = 1;
    while (N < (x << 1)) N <<= 1;
    Inv = pow_mod(N, MOD - 2);
    for (int i = 1; i < N; i++)
        if (i & 1)
            rev[i] = (rev[i >> 1] >> 1) | (N >> 1);
        else
            rev[i] = (rev[i >> 1] >> 1);
}
inline int Calc(int x)
{
    N = 1;
    while (N < (x << 1)) N <<= 1;
    return N;
}
void FFt(int *a, int op)
{
    int w, wn, t;
    for (int i = 0; i < N; i++)
        if (i < rev[i])
            swap(a[i], a[rev[i]]);
    for (int k = 2; k <= N; k <<= 1)
    {
        wn = pow_mod(3, op == 1 ? (MOD - 1) / k : MOD - 1 - (MOD - 1) / k);
        for (int j = 0; j < N; j += k)
        {
            w = 1;
            for (int i = 0; i < (k >> 1); i++, w = 1ll * w * wn % MOD)
            {
                t = 1ll * a[i + j + (k >> 1)] * w % MOD;
                a[i + j + (k >> 1)] = (a[i + j] - t + MOD) % MOD;
                a[i + j] = (a[i + j] + t) % MOD;
            }
        }
    }
    if (op == -1)
        for (int i = 0; i < N; i++)
            a[i] = 1ll * a[i] * Inv % MOD;
}
int tmp[MAXN], x[MAXN];
void Get_Inv(int *a, int *b, int n)
{
    if (n == 1) return b[0] = pow_mod(a[0], MOD - 2), void();
    Get_Inv(a, b, n + 1 >> 1);
    Init(n);
    for (int i = 0; i < n; i++) tmp[i] = a[i];
    for (int i = n; i < N; i++) tmp[i] = 0;
    FFt(tmp, 1), FFt(b, 1);
    for (int i = 0; i < N; i++) 
        b[i] = 1ll * b[i] * ((2 - 1ll * b[i] * tmp[i] % MOD + MOD) % MOD) % MOD;
    FFt(b, -1);
    for (int i = n; i < N; i++) b[i] = 0;
}
int Get_mod(int *a, int ra, int *b, int rb, int *c)
{
    while (ra && !a[ra - 1]) --ra;
    while (rb && !b[rb - 1]) --rb;
    if (ra < rb)
    {
        memcpy (c, a, ra << 2);
        memset (c + ra, 0, (rb - ra) << 2);
        return rb;
    }
    static int tmp1[MAXN], tmp2[MAXN];
    int rc = ra - rb + 1;
    int l = Calc(rc);
    for (int i = 0; i < l; i++) tmp1[i] = 0;
    reverse_copy(b, b + rb, tmp1);
    for (int i = 0; i < l; i++) tmp2[i] = 0; 
    Get_Inv(tmp1, tmp2, rc);
    for (int i = rc; i < l; i++) tmp2[i] = 0;
    reverse_copy(a, a + ra, tmp1);
    for (int i = rc; i < l; i++) tmp1[i] = 0;
    Init(rc), FFt(tmp1, 1), FFt(tmp2, 1);
    for (int i = 0; i < N; i++) tmp1[i] = 1ll * tmp1[i] * tmp2[i] % MOD;
    FFt(tmp1, -1); 
    reverse(tmp1, tmp1 + rc); 
    Init(ra);
    for (int i = 0; i < rb; i++) tmp2[i] = b[i];
    for (int i = rb; i < N; i++) tmp2[i] = 0;
    for (int i = rc; i < N; i++) tmp1[i] = 0;
    FFt(tmp1, 1), FFt(tmp2, 1);
    for (int i = 0; i < N; i++) tmp1[i] = 1ll * tmp1[i] * tmp2[i] % MOD;
    FFt(tmp1, -1);
    for (int i = 0; i < rb; i++) c[i] = (a[i] - tmp1[i] + MOD) % MOD;
    for (int i = rb; i < N; i++) c[i] = 0;
    while (rb && !c[rb - 1]) --rb;
    return rb;
}
int Solve0(int *a, int l, int r)
{
    int ra = r - l + 2;
    if (ra <= 256)
    {
        memset(a, 0, ra << 2), a[0] = 1;
        for (int i = l; i <= r; i++)
            for (int v = x[i], j = i - l; j != -1; j--)
            {
                a[j + 1] = (a[j + 1] + a[j]) % MOD;
                a[j] = 1ll * a[j] * v % MOD;
            }
        return ra;
    }
    int mid = l + r >> 1;
    int *f1 = a, r1 = Solve0(f1, l, mid);
    int *f2 = a + r1, r2 = Solve0(f2, mid + 1, r);
    N = 1;
    while (N < ra) N <<= 1;
    Inv = pow_mod(N, MOD - 2);
    for (int i = 1; i < N; i++)
        if (i & 1)
            rev[i] = (rev[i >> 1] >> 1) | (N >> 1);
        else
            rev[i] = (rev[i >> 1] >> 1);
    static int tmp1[L], tmp2[L];
    memcpy(tmp1, f1, r1 << 2), memset (tmp1 + r1, 0, (N - r1) << 2), FFt(tmp1, 1);
    memcpy(tmp2, f2, r2 << 2), memset (tmp2 + r2, 0, (N - r2) << 2), FFt(tmp2, 1);
    for (int i = 0; i < N; i++) a[i] = 1ll * tmp1[i] * tmp2[i] % MOD;
    FFt(a, -1);
    return ra;
}
int *p[L];
int sta[MAXN];
int mem[LM << 4], *head = mem;
inline int Solve1(int id, int l, int r)
{
    int ra = r - l + 2;
    if (ra <= 256)
    {
        int *f = p[id] = head; head += ra;
        memset (f, 0, ra << 2), 0[f] = 1;
        for (int i = l; i <= r; i++)
            for (int v = (MOD - sta[i]) % MOD, j = i - l; j != -1; j--)
                f[j + 1] = (f[j + 1] + f[j]) % MOD, f[j] = 1ll * f[j] * v % MOD;
        return ra;
    }
    int mid = l + r >> 1;
    int r1 = Solve1(id << 1, l, mid), *f1 = p[id << 1];
    int r2 = Solve1(id << 1 | 1, mid + 1, r), *f2 = p[id << 1 | 1];
    N = 1;
    while (N < ra) N <<= 1;
    Inv = pow_mod(N, MOD - 2);
    for (int i = 1; i < N; i++)
        if (i & 1)
            rev[i] = (rev[i >> 1] >> 1) | (N >> 1);
        else
            rev[i] = (rev[i >> 1] >> 1);
    static int tmp1[LM], tmp2[LM];
    memcpy(tmp1, f1, r1 << 2), memset (tmp1 + r1, 0, (N - r1) << 2), FFt(tmp1, 1);
    memcpy(tmp2, f2, r2 << 2), memset (tmp2 + r2, 0, (N - r2) << 2), FFt(tmp2, 1);
    int *f = p[id] = head; head += ra;
    for (int i = 0; i < N; i++) f[i] = 1ll * tmp1[i] * tmp2[i] % MOD;
    FFt(f, -1);
    return ra; 
}
int val0[LM], val[LM];
void Solve2(int id, int *a, int *b, int l, int r, int deg)
{
    int ra = r - l + 2;
    if (deg <= 256)
    {
        int F, G;
        for (int i = l; i <= r; i++)
        {
            F = G = 0;
            int u = sta[i], v = 1;
            for (int j = 0; j <= deg; j++)
            {
                F = (F + 1ll * v * a[j]) % MOD;
                G = (G + 1ll * v * b[j]) % MOD;
                v = 1ll * v * u % MOD;
            }
            val0[i] = F, val[i] = G;
        }
        return;
    }
    int mid = l + r >> 1;
    int r1 = Get_mod(a, deg, p[id], ra, a + deg); a += deg;
    int r2 = Get_mod(b, deg, p[id], ra, b + deg); b += deg;
    ra = min(r1, r2);
    Solve2(id << 1, a, b, l, mid, ra), Solve2(id << 1 | 1, a, b, mid + 1, r, ra);
}
int g[MAXN], h[MAXN];
int main()
{
    n = read(), m = read();
    Sum = 0;
    for (int i = 1; i <= n; i++)
        x[i] = read() % MOD, Sum = (Sum + x[i]) % MOD;
    int A = 1, B = 0, C = 0, D = 1;
    for (int i = 1; i <= m; i++)
    {
        if (read() == 1)
        {
            int v = read() % MOD;
            A = (A + 1ll * v * B % MOD) % MOD;
            C = (C + 1ll * v * D % MOD) % MOD;
            insert(A, B, C, D, i);
        }
        else
        {
            swap(A, B);
            swap(C, D);
            insert(A, B, C, D, i);
        }
    }
    if (cnt)
    {
        for (int i = 1; i <= cnt; i++) sta[i] = s[i].g;
        sort(sta + 1, sta + cnt + 1);
        int lenth = unique(sta + 1, sta + cnt + 1) - sta - 1;
        for (int i = 1; i <= cnt; i++)
            s[i].g = lower_bound(sta + 1, sta + lenth + 1, s[i].g) - sta;
        Solve0(g, 1, n);
        for (int i = 1; i <= n; i++) h[i - 1] = 1ll * g[i] * i % MOD;
        Solve1(1, 1, lenth);
        Solve2(1, g, h, 1, lenth, n + 1);
        for (int i = 1; i <= cnt; i++)
            ans[s[i].id] = (1ll * s[i].e * n % MOD + 1ll * s[i].f * val[s[i].g] % MOD * pow_mod(val0[s[i].g], MOD - 2) % MOD) % MOD;
    }
    for (int i = 1; i <= m; i++)
        printf ("%d\n", ans[i]);
    // while (1);
}
```

##### 高斯消元
```c++
long long gauss(int n)
{
    long long ans = 1;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            a[i][j] = (a[i][j] + MOD) % MOD;
    for (int i = 1; i <= n; i++)
    {
        for (int j = i + 1; j <= n; j++)
            while (a[j][i])
            {
                int t = a[i][i] / a[j][i];
                for (int k = i; k <= n; k++)
                {
                    a[i][k] = (a[i][k] - a[j][k] * t % MOD + MOD) % MOD;
                    swap(a[i][k], a[j][k]);
                }
                ans = (MOD - ans) % MOD;
            }
        ans = ans * a[i][i] % MOD;
    }
    return ans;
}



void gauss() {
  int im, num = 1;
  for (int k = 1; k <= n; k++, num++) {
    im = k;
    for (int i = k + 1; i <= n; i++) {
      if (fabs(a[i][k]) > fabs(a[im][k])) i = im;
    }
    if (im != k) {
      for (int i = k; i <= n + 1; i++) swap(a[num][i], a[im][i]);
    }
    if (!a[num][k]) {
      num--;
      continue;
    }
    for (int i = num + 1; i <= n; i++) {
      if (!a[num][k]) continue;
      long double t = a[i][k] / a[num][k];
      for (int j = k; j <= n + 1; j++) {
        a[i][j] -= t * a[k][j];
      }
    }
  }
  for (int i = n; i >= 1; i--) {
    for (int j = n; j >= i + 1; j--) {
      a[i][n + 1] -= a[i][j] * x[j];
    }
    x[i] = a[i][n + 1] / a[i][i];
    sum += x[i];
  }
}
```

#### BSGS&exBSGS

```c++
int pow_mod(int a, int b, int c)
{
    int ans = 1;
    while (b)
    {
        if (b & 1) ans = ans * a % c;
        b >>= 1;
        a = a * a % c;
    }
    return ans;
}
struct Hash_Table
{
    struct edge
    {
        int next, x, ans;
    }v[100005];
    int first[76545], p;
    int &operator[](int x)
    {
        int H = x % 76543;
        for (int i = first[H]; i != -1; i = v[i].next)
            if (v[i].x == x)
                return v[i].ans;
        v[p].x = x;
        v[p].next = first[H];
        first[H] = p++;
        return v[p - 1].ans = 0;
    }
    bool count(int x)
    {
        int H = x % 76543;
        for (int i = first[H]; i != -1; i = v[i].next)
            if (v[i].x == x)
                return 1;
        return 0;
    }
    void clear()
    {
        memset (first, -1, sizeof (first));
        p = 0;
    }
}Hash;
int gcd(int a, int b)
{
    return b == 0 ? a : gcd(b, a % b);
}
int BSGS(int a, int b, int c, int d)
{
    Hash.clear();
    int m = floor(sqrt(c)) + 1;
    for (int i = 0; i <= m; i++)
    {
        // if (!Hash.count(b)) 
        Hash[b] = i;
        b = b * a % c;
    }
    int s = pow_mod(a, m, c);
    for (int i = 1; i <= m; i++)
    {
        d = d * s % c;
        if (Hash.count(d)) return i * m - Hash[d];
    }
    return -1;
}
int exBSGS(int a, int b, int c)
{
    b %= c;
    int s = 1;
    for (int i = !a; i <= 40; i++)
    {
        if (s == b) return i;
        s = s * a % c;
    }
    int cnt = 0, d = 1;
    for (int i; (i = gcd(a, c)) != 1; cnt++)
    {
        if (b % i) return -1;
        b /= i, c /= i;
        d = d * a / i % c;
    }
    int ret = BSGS(a, b, c, d);
    if (ret == -1) return -1;
    return ret + cnt;
}
```