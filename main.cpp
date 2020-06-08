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
/*#include<stdio.h>
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
   int W;
   int n = sizeof(val)/sizeof(val[0]);


   cout << "Enter capacity of knapsack: ";
   cin >> W;

   cout << "The solution is : " << knapsack(W, wt, val, n);
   return 0;
}include<stdio.h>
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
   int W;
   int n = sizeof(val)/sizeof(val[0]);
include<stdio.h>
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
   int W;
   int n = sizeof(val)/sizeof(val[0]);


   cout << "Enter capacity of knapsack: ";
   cin >> W;

   cout << "The solution is : " << knapsack(W, wt, val, n);
   return 0;
}

   cout << "Enter capacity of knapsack: ";
   cin >> W;

   cout << "The solution is : " << knapsack(W, wt, val, n);
   return 0;
}*/

/*#include <iostream>
using namespace std;
int max(int x, int y) {
   return (x > y) ? x : y;
}
int knapSack(int W, int w[], int v[], int n) {
   int i, wt;
   int K[n + 1][W + 1];
   for (i = 0; i <= n; i++) {
      for (wt = 0; wt <= W; wt++) {
         if (i == 0 || wt == 0)
         K[i][wt] = 0;
         else if (w[i - 1] <= wt)
            K[i][wt] = max(v[i - 1] + K[i - 1][wt - w[i - 1]], K[i - 1][wt]);
         else
        K[i][wt] = K[i - 1][wt];
      }
   }
   return K[n][W];
}
int main() {
   cout << "Enter the number of items in a Knapsack: ";
   int n, W;
   cin >> n;
   int v[n], w[n];
   for (int i = 0; i < n; i++) {
      cout << "\nEnter value for item " << i << ": ";
      cin >> v[i];
      cout << "Enter weight for item " << i << ": ";
      cin >> w[i];
   }
   cout << "\nEnter the capacity of knapsack: ";
   cin >> W;
   cout << "\nThe solution is : "  <<knapSack(W, w, v, n);
   return 0;
}*/

class Simplex{

    private:
        int rows, cols;
        //stores coefficients of all the variables
        std::vector <std::vector<float> > A;
        //stores constants of constraints
        std::vector<float> B;
        //stores the coefficients of the objective function
        std::vector<float> C;

        float maximum;

        bool isUnbounded;

    public:
        Simplex(std::vector <std::vector<float> > matrix,std::vector<float> b ,std::vector<float> c){
            maximum = 0;
            isUnbounded = false;
            rows = matrix.size();
            cols = matrix[0].size();
            A.resize( rows , vector<float>( cols , 0 ) );
            B.resize(b.size());
            C.resize(c.size());




            for(int i= 0;i<rows;i++){             //pass A[][] values to the metrix
                for(int j= 0; j< cols;j++ ){
                    A[i][j] = matrix[i][j];

                }
            }



            for(int i=0; i< c.size() ;i++ ){      //pass c[] values to the B vector
                C[i] = c[i] ;
            }
            for(int i=0; i< b.size();i++ ){      //pass b[] values to the B vector
                B[i] = b[i];
            }




        }

        bool simplexAlgorithmCalculataion(){
            //check whether the table is optimal,if optimal no need to process further
            if(checkOptimality()==true){
                            return true;
            }

            //find the column which has the pivot.The least coefficient of the objective function(C array).
            int pivotColumn = findPivotColumn();


            if(isUnbounded == true){
                cout<<"Error unbounded"<<endl;
                            return true;
            }

            //find the row with the pivot value.The least value item's row in the B array
            int pivotRow = findPivotRow(pivotColumn);

            //form the next table according to the pivot value
            doPivotting(pivotRow,pivotColumn);

            return false;
        }

        bool checkOptimality(){
             //if the table has further negative constraints,then it is not optimal
            bool isOptimal = false;
            int positveValueCount = 0;

            //check if the coefficients of the objective function are negative
            for(int i=0; i<C.size();i++){
                float value = C[i];
                if(value >= 0){
                    positveValueCount++;
                }
            }
            //if all the constraints are positive now,the table is optimal
            if(positveValueCount == C.size()){
                isOptimal = true;
                print();
            }
            return isOptimal;
        }

        void doPivotting(int pivotRow, int pivotColumn){

            float pivetValue = A[pivotRow][pivotColumn];//gets the pivot value

            float pivotRowVals[cols];//the column with the pivot

            float pivotColVals[rows];//the row with the pivot

            float rowNew[cols];//the row after processing the pivot value

            maximum = maximum - (C[pivotColumn]*(B[pivotRow]/pivetValue));  //set the maximum step by step
             //get the row that has the pivot value
             for(int i=0;i<cols;i++){
                pivotRowVals[i] = A[pivotRow][i];
             }
             //get the column that has the pivot value
             for(int j=0;j<rows;j++){
                pivotColVals[j] = A[j][pivotColumn];
            }

            //set the row values that has the pivot value divided by the pivot value and put into new row
             for(int k=0;k<cols;k++){
                rowNew[k] = pivotRowVals[k]/pivetValue;
             }

            B[pivotRow] = B[pivotRow]/pivetValue;


             //process the other coefficients in the A array by subtracting
             for(int m=0;m<rows;m++){
                //ignore the pivot row as we already calculated that
                if(m !=pivotRow){
                    for(int p=0;p<cols;p++){
                        float multiplyValue = pivotColVals[m];
                        A[m][p] = A[m][p] - (multiplyValue*rowNew[p]);
                        //C[p] = C[p] - (multiplyValue*C[pivotRow]);
                        //B[i] = B[i] - (multiplyValue*B[pivotRow]);
                    }

                }
             }

            //process the values of the B array
            for(int i=0;i<B.size();i++){
                if(i != pivotRow){

                        float multiplyValue = pivotColVals[i];
                        B[i] = B[i] - (multiplyValue*B[pivotRow]);

                }
            }
                //the least coefficient of the constraints of the objective function
                float multiplyValue = C[pivotColumn];
                //process the C array
                 for(int i=0;i<C.size();i++){
                    C[i] = C[i] - (multiplyValue*rowNew[i]);

                }
             //replacing the pivot row in the new calculated A array
             for(int i=0;i<cols;i++){
                A[pivotRow][i] = rowNew[i];
             }
        }

        //print the current A array
        void print(){
            for(int i=0; i<rows;i++){
                for(int j=0;j<cols;j++){
                    cout<<A[i][j] <<" ";
                }
                cout<<""<<endl;
            }
            cout<<""<<endl;
        }

        //find the least coefficients of constraints in the objective function's position
        int findPivotColumn(){

            int location = 0;
            float minm = C[0];



            for(int i=1;i<C.size();i++){
                if(C[i]<minm){
                    minm = C[i];
                    location = i;
                }
            }

            return location;

        }

        //find the row with the pivot value.The least value item's row in the B array
        int findPivotRow(int pivotColumn){
            float positiveValues[rows];
            std::vector<float> result(rows,0);
            //float result[rows];
            int negativeValueCount = 0;

            for(int i=0;i<rows;i++){
                if(A[i][pivotColumn]>0){
                    positiveValues[i] = A[i][pivotColumn];
                }
                else{
                    positiveValues[i]=0;
                    negativeValueCount+=1;
                }
            }
            //checking the unbound condition if all the values are negative ones
            if(negativeValueCount==rows){
                isUnbounded = true;
            }
            else{
                for(int i=0;i<rows;i++){
                    float value = positiveValues[i];
                    if(value>0){
                        result[i] = B[i]/value;

                    }
                    else{
                        result[i] = 0;
                    }
                }
            }
            //find the minimum's location of the smallest item of the B array
            float minimum = 99999999;
            int location = 0;
            for(int i=0;i<sizeof(result)/sizeof(result[0]);i++){
                if(result[i]>0){
                    if(result[i]<minimum){
                        minimum = result[i];

                        location = i;
                    }
                }

            }

            return location;

        }

        void CalculateSimplex()
        {
            bool end = false;

            cout<<"Initial matrix (Not the optimal solution)"<<endl;
            print();

            cout<<" "<<endl;
            cout<<"Final matrix (Optimal solution)"<<endl;


            while(!end)
            {
                bool result = simplexAlgorithmCalculataion();

                if(result==true)
                    {
                    end = true;
                    }
            }

            cout<<"Answers for the Constraints of variables"<<endl;

            for(int i=0;i< A.size(); i++){  //every basic column has the values, get it form B array
                int count0 = 0;
                int index = 0;
                for(int j=0; j< rows; j++){
                    if(A[j][i]==0.0){
                        count0 += 1;
                    }
                    else if(A[j][i]==1){
                        index = j;
                    }
                }

                if(count0 == rows -1 )
                {

                    cout<<"variable"<<index+1<<": "<<B[index]<<endl;  //every basic column has the values, get it form B array
                }

                else
                {
                    cout<<"variable"<<index+1<<": "<<0<<endl;

                }
            }
           cout<<"  "<<endl;
           cout<<"maximum value: "<<maximum<<endl;  //print the maximum values
        }
};

int main()
{
    int colSizeA;
    cout << "Column size (slack variables included):  ";
    cin >> colSizeA;   //should initialise columns size in A

    int rowSizeA;  //should initialise columns row in A[][] vector
    cout << "Row size: ";
    cin >> rowSizeA;

   // float C[]= {-6,-5,-4,0,0,0};  //should initialis the c arry here
    float B[]={5000,5000,4800};  // should initialis the b array here

    int objective;
    cout << "How many variables: ";
    cin >> objective;

    float C[objective];
    cout << "Enter the objective function: ";
    for(int f=0; f<objective; f++)
    {
          cin >> C[f];
    }

    float a[rowSizeA][colSizeA];

    cout << "Enter the elements of the matrix: " <<endl;
    for(int w=0; w<rowSizeA; w++)
    {
        for(int q=0; q<colSizeA; q++)
        {
          cin >> a[w][q];
        }
    }
 /*   for(int s=1; s<=rowSizeA; s++)
    {
        for (int d=1; d<=colSizeA; d++)
        {
            if (s==d)
            {
                cout << 1;
            }
            else
            {
             cout << 0;
            }
        }
        cout << endl;
    }*/

        std::vector <std::vector<float> > vec2D(rowSizeA, std::vector<float>(colSizeA, 0));

        std::vector<float> b(rowSizeA,0);
        std::vector<float> c(colSizeA,0);




       for(int i=0;i<rowSizeA;i++){         //make a vector from given array
            for(int j=0; j<colSizeA;j++){
                vec2D[i][j] = a[i][j];
            }
       }





       for(int i=0;i<rowSizeA;i++){
            b[i] = B[i];
       }

        for(int i=0;i<colSizeA;i++){
            c[i] = C[i];
       }


      // hear the make the class parameters with A[m][n] vector b[] vector and c[] vector
      Simplex simplex(vec2D,b,c);
      simplex.CalculateSimplex();


    return 0;
}
