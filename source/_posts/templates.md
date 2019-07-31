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
	1. 树状数组 (BIT)
	1. 线段树 (SGM)
		1. 划分树与归并树
	1. 并查集 (UFS)
		1. 带权并查集
		1. 路径压缩
		1. 按秩合并
	1. Sparse Table (ST)
	1. K维树 (KDT)
	1. 动态树
		1. 点/边分治树 (DCT)
		1. Link-Cut Tree (LCT)
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
		1. dijkstra 算法 (DJ)
		1. Bellman-Ford 算法 (BFD)
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
			1. 高斯消元
			1. 拟阵
			1. Matrix-Tree 定理
			1. 线性递推
	1. 多项式与幂级数
		1. [DFT/FFT (FT)](#FFT)
			1. [NTT](#NTT)
			1. Bluestein 算法
		1. 多项式基本运算
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
			1. BSGS
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