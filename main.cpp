#include <iostream>
#include <limits>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <queue>
#include <bits/stdc++.h>
/*using namespace std;

#define TABSIZE 4

void FayardPlateau(float* content, float* cost, float* utility, float maxCost, unsigned int size);
void BranchAndBound(bool* content, float* cost, float* utility, float maxCost, unsigned int size);


 *

int main(){
	bool content[TABSIZE] = {
		false,
		false,
		false,
		false
	};

	float utility[TABSIZE]	=	{8.f, 18.f, 20.f, 11.f};
	float cost[TABSIZE]		=	{3.f, 7.f, 9.f, 6.f};

	BranchAndBound(content, cost, utility, 17.f, TABSIZE);

	for(int i = 0; i < TABSIZE; i++){
		std::cout << content[i] << endl;
	}std::cout << '\n';

	float fcost = 0.f;
	for(int i = 0; i < TABSIZE; i++){
		fcost += content[i] * cost[i];
	}std::cout << "Cost = " << fcost << endl;

	return 0;
}



struct uc{
	float ucRatio;
	unsigned int index;
};

int cmp(const void* elm1, const void* elm2){
	int a = ((uc*)elm1)->ucRatio;
	int b = ((uc*)elm2)->ucRatio;
	return a > b ? -1 : a < b;
}


 //*	BRANCH & BOUND


struct bbNode{
	unsigned int	lev;
	bool*			set;
	float			benef;
	float			cost;
	float			bound;
};

float calculateBound(bbNode nd, float* cost, float* utility, uc* ucTab, float maxCost, unsigned int size){
	float addCost = 0.f;
	float addBenef = 0.f;

	for(int i = nd.lev; i < size; i++){
		if(nd.cost + addCost + cost[ucTab[i].index] > maxCost){
			//	Partial add of the best c/u object
			float partialBenefit = ((maxCost - addCost) / cost[ucTab[i].index]) * utility[ucTab[i].index];
			return nd.benef + addBenef + partialBenefit;
		}
		addCost += cost[ucTab[i].index];
		addBenef += utility[ucTab[i].index];
	}
	return nd.cost + addCost;
}

float calculateBenefit(bbNode nd, float* utility, unsigned int size){
	float res = 0.f;
	for(int i = 0; i < size; i++)
		if(nd.set[i])
			res+=utility[i];
	return res;
}

void cpyNode(bbNode* dest, bbNode src, unsigned int size){
	*dest = src;
	dest->set = (bool*)malloc(size * sizeof(bool));
	memcpy(dest->set, src.set, size * sizeof(bool));
}

void setNextObject(bbNode* nd, bool val, float* cost, uc* ucTab){
	nd->set[ucTab[nd->lev].index] = val;
	nd->cost += val * cost[ucTab[nd->lev].index];
	nd->lev++;
}

void BranchAndBound(bool* content, float* cost, float* utility, float maxCost, unsigned int size){
	//	Declarations
	uc* ucTab = (uc*)malloc(size * sizeof(uc));
	queue<bbNode> queueNd;

	//	Sorting items by Cost/Utility ratio
	for(int i = 0; i < size; i++){
		ucTab[i].ucRatio = utility[i] / cost[i];
		ucTab[i].index = i;
	}qsort(ucTab, size, sizeof(uc), cmp);

	//	Initializing level, set of items selected etc.
	bbNode bsf;

	bsf.lev = 0;
	bsf.set = (bool*)malloc(size * sizeof(bool));
	bsf.benef = - INFINITY;
	bsf.cost = 0.f;
	bsf.bound = 0.f;
	memset(bsf.set, 0, size * sizeof(bool));

	bbNode aNode;

	aNode.lev = 0;
	aNode.set = (bool*)malloc(size * sizeof(bool));
	aNode.benef = 0.f;
	aNode.cost = 0.f;
	aNode.bound = calculateBound(bsf, cost, utility, ucTab, maxCost, size);
	memset(aNode.set, 0, size * sizeof(bool));

	//	Enqueuing
	queueNd.push(aNode);

	while(!queueNd.empty()){
		bbNode currNode = queueNd.front();
		queueNd.pop();

		memcpy(content, currNode.set, size * sizeof(bool));

		float _cost = 0.f;
		cout << "Contenu : \n";
		for(int i = 0; i < size; i++){
			_cost += currNode.set[i] * cost[i];
			cout << currNode.set[i] << endl;
		}cout << "Cout : " << _cost << endl << endl;

		if(currNode.bound > bsf.benef){
			//	Create a node nextAdded equal to currNode with the next item added
			bbNode nextAdded;
			cpyNode(&nextAdded, currNode, size);
			delete[] currNode.set;	//	Memory management
			setNextObject(&nextAdded, true, cost, ucTab);
			nextAdded.bound = calculateBound(nextAdded, cost, utility, ucTab, maxCost, size);
			nextAdded.benef = calculateBenefit(nextAdded, utility, size);

			//	if nextAdded's benefit > bestSoFar's benefit
			if(nextAdded.benef > bsf.benef){
				//	set bsf equal to nextAdded
				delete[] bsf.set;
				cpyNode(&bsf, nextAdded, size);
			}

			//	if nextAdded's bound > bestSoFar's benefit
			if(nextAdded.bound > bsf.benef){
				//	enqueue nextAdded
				queueNd.push(nextAdded);
			}
		}

		//	Create a node nextNotAdded equal to currNode without the next item added
		bbNode nextNotAdded;
		cpyNode(&nextNotAdded, currNode, size);
		delete[] currNode.set;	//	Memory management
		setNextObject(&nextNotAdded, false, cost, ucTab);

		nextNotAdded.bound = calculateBound(nextNotAdded, cost, utility, ucTab, maxCost, size);
		nextNotAdded.benef = calculateBenefit(nextNotAdded, utility, size);

		//	if nextNotAdded's bound > bsf's benef
		if(nextNotAdded.bound > bsf.benef){
			//	enqueue nextNotAdded
			queueNd.push(nextNotAdded);
		}
	}
}



 //*	FAYARD PLATEAU



void FayardPlateau(float* content, float* cost, float* utility, float maxCost, unsigned int size){
	//	Tri des couts/utilité
	uc* ucTab = new uc[size];

	for(int i = 0; i < size; i++){
		ucTab[i].ucRatio = utility[i] / cost[i];
		ucTab[i].index = i;
	}

	qsort(ucTab, size, sizeof(uc), cmp);

	//	Ajout des objets dans le sac
	float totalCost = 0.f;
	for(int i = 0; i < size; i++){
		if(cost[ucTab[i].index] <= maxCost - totalCost){	//	Ajout total
			content[ucTab[i].index] = 1.f;
			totalCost += cost[ucTab[i].index];
		}
		else{												//	Ajout partiel
			float ratioAdd = (maxCost - totalCost) / cost[ucTab[i].index];
			content[ucTab[i].index] = ratioAdd;
			totalCost += cost[ucTab[i].index] * ratioAdd;
		}
	}
}*/
/*sing namespace std;

struct Item
{
    float weight;
    int value;
};
struct Node
{
    int level, profit, bound;
    float weight;
};

bool cmp(Item a, Item b)
{
    double r1 = (double)a.value / a.weight;
    double r2 = (double)b.value / b.weight;
    return r1 > r2;
}
int bound(Node u, int n, int W, Item arr[])
{
    if (u.weight >= W)
        return 0;
    int profit_bound = u.profit;
    int j = u.level + 1;
    int totweight = u.weight;

    while ((j < n) && (totweight + arr[j].weight <= W))
    {
        totweight    = totweight + arr[j].weight;
        profit_bound = profit_bound + arr[j].value;
        j++;
    }
    if (j < n)
        profit_bound = profit_bound + (W - totweight) * arr[j].value /
                                         arr[j].weight;

    return profit_bound;
}

int knapsack(int W, Item arr[], int n)
{
    sort(arr, arr + n, cmp);
    queue<Node> Q;
    Node u, v;
    u.level = -1;
    u.profit = u.weight = 0;
    Q.push(u);
    int maxProfit = 0;
    while (!Q.empty())
    {
        u = Q.front();
        Q.pop();
        if (u.level == -1)
            v.level = 0;

        if (u.level == n-1)
            continue;
        v.level = u.level + 1;
        v.weight = u.weight + arr[v.level].weight;
        v.profit = u.profit + arr[v.level].value;
        if (v.weight <= W && v.profit > maxProfit)
            maxProfit = v.profit;
        v.bound = bound(v, n, W, arr);
        if (v.bound > maxProfit)
            Q.push(v);
        v.weight = u.weight;
        v.profit = u.profit;
        v.bound = bound(v, n, W, arr);
        if (v.bound > maxProfit)
            Q.push(v);
    }

    return maxProfit;
}
int main()
{
    int W = 28;   // Weight of knapsack
    Item arr[] = {{1, 1}, {4, 7}};
    int n = sizeof(arr) / sizeof(arr[0]);

    cout << "Maximum possible profit = "
         << knapsack(W, arr, n);

    return 0;
}*/
#include <iostream>
#include <bits/stdc++.h>
using namespace std;
/*typedef struct
{
   int v;
   int w;
   float d;
} Item;

void input(Item items[],int sizeOfItems)
{
   cout<< "Enter size: ";
   cin >> sizeOfItems;
   cout << "Enter total "<< sizeOfItems <<" item's values and weight" <<
   endl;
   for(int i = 0; i < sizeOfItems; i++)
   {
      cout << "Enter "<< i+1 << " V: ";
      cin >> items[i].v;
      cout << "Enter "<< i+1 << " W: ";
      cin >> items[i].w;
   }
}

void display(Item items[], int sizeOfItems)
{
   int i;
   cout << "values: ";
   for(i = 0; i < sizeOfItems; i++)
   {
      cout << items[i].v << "\t";
   }
   cout << endl << "weight: ";
   for (i = 0; i < sizeOfItems; i++)
   {
      cout << items[i].w << "\t";
   }
   cout << endl;
}

bool compare(Item i1, Item i2)
{
   return (i1.d > i2.d);
}

float knapsack(Item items[], int sizeOfItems, int W)
{
   int i, j, pos;
   Item mx, temp;
   float totalValue = 0, totalWeight = 0;
   for (i = 0; i < sizeOfItems; i++)
   {
      items[i].d = items[i].v / items[i].w;
   }
   sort(items, items+sizeOfItems, compare);
   for(i=0; i<sizeOfItems; i++)
    {
      if(totalWeight + items[i].w<= W)
      {
         totalValue += items[i].v ;
         totalWeight += items[i].w;
      }
   else
      {
         int wt = W-totalWeight;
         totalValue += (wt * items[i].d);
         totalWeight += wt;
         break;
      }
   }
   cout << "total weight in bag " << totalWeight<<endl;
   return totalValue;
}

int main()
{
   int W;
   Item items[4];
   input(items, 4);
   cout << "Entered data \n";
   display(items,4);
   cout<< "Enter Knapsack weight \n";
   cin >> W;
   float mxVal = knapsack(items, 4, W);
   cout << "Max value for "<< W <<" weight is "<< mxVal;
}*/
#include<stdio.h>
int max(int a, int b) {
   if(a>b){
      return a;
   } else {
      return b;
   }
}
int knapsack(int W, int wt[], int val[], int n) {
   int i, w;
   int knap[n+1][W+1];
   for (i = 0; i <= n; i++) {
      for (w = 0; w <= W; w++) {
         if (i==0 || w==0)
            knap[i][w] = 0;
         else if (wt[i-1] <= w)
            knap[i][w] = max(val[i-1] + knap[i-1][w-wt[i-1]], knap[i-1][w]);
         else
            knap[i][w] = knap[i-1][w];
      }
   }
   return knap[n][W];
}
int main() {
   int val[] = {40, 30, 50, 10};
   int wt[] = {2, 5, 10, 5};
   int W = 16;
   int n = sizeof(val)/sizeof(val[0]);
   cout << "The solution is : " << knapsack(W, wt, val, n);
   return 0;
}
