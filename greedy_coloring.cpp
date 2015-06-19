///input===> pairs.txt; omit the first 6 lines; from 0 to N-1;///////////////
///output===> Nb vs. lb(lb*)

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>

using namespace std;

const int MAX=2501;
const int MAXL=10;
const int MINL=3;
int N=0;
struct network
{
    int b;//unburnt--->0; burnt--->1
    vector<int> neighbor;
    int geodis[MAX];
    int boxid;  // =color value; start from 1
};
network node[MAX];

void readnet(network *);
void printnet(network *, int);
void burning(network *);//calculate node[i].geodis[j]; return averaged path

bool complete(network *); // determine whether all nodes are covered
float boxcr(network *,int);//return lB* --- averaged real lB
float module(network *,int);//calculate Modularity
int GC(network * node, int lb);//greedy coloring

int main()
{
    srand(time(NULL));   
    
    ofstream outFile;
    outFile.open("GC.txt");

    outFile << "#lb" << "\t" << "lb*" << "\t" << "M" << "\t" << "Nb" << endl;

    
    const char Filename[1][40]=
{
    "partition_l9.clu",
    //"partition_l2.clu",
    //"partition_l3.clu",
    //"partition_l4.clu",
    //"partition_l5.clu",
};

    readnet(node); // N=# of nodes
    
    burning(node);
    cout << "burning Done!" << endl;

    int Nb;
    float lb_cr;
    float Mb;
    for (int lb=MINL; lb<=MAXL; lb++)
    {
        Nb=GC(node,lb);
        lb_cr=boxcr(node,Nb);// lB*
        Mb=module(node,Nb);// Modularity
        
        if (lb==9)
        {
            ofstream pajek;
            pajek.open(Filename[lb-MINL]);
            pajek << "*vertices " << N << endl;
            for (int i=0;i<N;i++)
                pajek << node[i].boxid << endl;
            pajek.close();
        }
        
        outFile << lb << "\t" << lb_cr << "\t" << Mb << "\t" << Nb << endl;
    }
    outFile.close();
    return 0;
}
// Greedy Coloring ---> cal node[i].boxid & return Nb for a given lb

int GC(network * node, int lb)
{
    node[0].boxid=1;
    for (int i=1;i<N;i++)
        node[i].boxid = 0; // mark all nodes(but node[0]) as uncolored

    int Nb=1;//---current color!!!
    int p;
    int r;
    vector<int> backup;
    for (int i=1 ; i<N ; i++)
    {
        backup.clear();//---backup color
        for (unsigned int q=0; q<node[i].neighbor.size(); q++)
        {
            p = node[i].neighbor[q];
            if (node[p].boxid > 0)
                backup.push_back(node[p].boxid);
        }
        for (int j=0 ; j<i ; j++)
            if (node[i].geodis[j] >= lb){
                for (unsigned int x=0; x<backup.size(); x++)
                    if (backup[x]==node[j].boxid)
                        backup.erase(backup.begin()+x);
            }
        /*
           cout << "node["<<i<<"]:";
           for (unsigned int t=0; t<opcl.size(); t++)
           cout << opcl[t] << ",";
           cout << endl;
        */
        if (backup.empty()){
            Nb++;
            node[i].boxid = Nb;
        }
        else {
            r = rand()%backup.size();
            node[i].boxid = backup[r];
        }
        
    }
    return Nb;
}

// burning ---> calculate node[i].geodis[j]

void burning(network * node)
{
    int i,t,c,g,h;
    unsigned int j,e,f,d;
    vector<int> active;
    vector<int> active0;
    bool flag;
    for (i=0;i<N;i++)
    {
        active.clear();
        active0.clear();
        for (c=0;c<N;c++)
        {
            node[c].b=0;//mark all nodes as unburnt
        }
        active.push_back(i);
        node[i].b=1;
        g=1;
        while (active.size()) 
        {
            active0=active;
            active.clear();
            for (j=0;j<active0.size();j++)
            {
                t=active0[j];
                for (d=0;d<node[t].neighbor.size();d++)
                {
                    h=node[t].neighbor[d];
                    if (!node[h].b)
                    {
                        flag=true;
                        for (e=0;e<active.size();e++)
                            if (active[e]==h)
                            {
                                flag=false;
                                break;
                            }
                        if(flag)
                        active.push_back(h);
                    }
                }
            }
            for (f=0;f<active.size();f++)
            {
                node[active[f]].b=1;
                node[i].geodis[active[f]]=g;
            }
            g++;
            //if (g>MAXL) break;-----no break in GC!
        }

    }

    /*
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++)
            cout << node[i].geodis[j] << " ";
    cout<<endl;}
    */
}

//---calculate Modularity

float module(network * node, int Nb)
{
    float out[Nb];
    float in[Nb];
    for (int i=0;i<Nb;i++)
    {
        out[i]=0;
        in[i]=0;
    }
    for (int i=0;i<N;i++)
        for (int j=i+1;j<N;j++)
            if (node[i].geodis[j]==1)
                if (node[i].boxid==node[j].boxid)
                    in[node[i].boxid-1]++;
                else 
                {
                     out[node[i].boxid-1]++;
                     out[node[j].boxid-1]++;
                }
    float m=0;
    for (int i=0;i<Nb;i++)
        m=m+in[i]/out[i];
    m=m/Nb;
    cout << m << endl;
    return m;
}


 //---box correction

float boxcr(network * node, int Nb)
{
    int lb[Nb];
    for (int t=0;t<Nb;t++)
    {
        int box[N];
        int mass=0;
        for (int u=0;u<N;u++)
            if (node[u].boxid==(t+1))
            {
                box[mass]=u;
                mass++;
            }
        //for (int i=0;i<mass;i++)
        //    cout << box[i];
        //cout << endl;
        lb[t]=0;
        for (int i=0;i<mass;i++)
            for (int j=i+1;j<mass;j++)
                if (node[box[i]].geodis[box[j]]>lb[t])
                    lb[t]=node[box[i]].geodis[box[j]];
    }
    //for (int i=0;i<Nb;i++)
    //    cout << lb[i] << endl;

    float lB=0;
    for (int i=0;i<Nb;i++)
        lB=lB+lb[i];
    lB=lB/Nb+1;
    //cout << lB << endl;
    return lB;
}

// --- read network & get N

void readnet(network * node)    
{
    //char Filename[] = "network.txt";
    char Filename[40];
    ifstream inFile;
    cout << "Please enter name of data file: ";
    cin.getline(Filename, 40);
    inFile.open(Filename);
    if (!inFile.is_open())
    {
        cout << "Could not open the file!" << Filename << endl;
        exit(EXIT_FAILURE);
    }
    for (int i=0;i<MAX;i++)
    {
        node[i].b = 0;// mark all nodes as unburnt --->0
        for (int j=0;j<MAX;j++)
        {
            node[i].geodis[j]=0;
        }  
    }
    string s;
    istringstream iss;

    for (int i=0; i<0; i++)//---!
    getline(inFile,s);

    while (getline(inFile,s))

    {
        int i,j;
        iss.str(s);
        iss>>i>>j;
        //i--;//when read .net file
        //j--;
        //cout << i << "\t<----->  \t" << j << endl;
        node[i].neighbor.push_back(j);
        node[j].neighbor.push_back(i);
        iss.clear();
    }
    for (int i=0;i<MAX;i++)
        if (node[i].neighbor.size())
            N++;
    cout <<  "# of nodes: " << N << endl;
    inFile.close();
}

// --- print out network

void printnet(network * node)    
{
    for (int i=0;i<N;i++) 
    {
        cout << "Node " << setw(2) << i << " has " << (int)node[i].neighbor.size(); 
        cout << " edges:";
        for (unsigned int j=0;j<node[i].neighbor.size();j++)
            cout << setw(3) << node[i].neighbor[j] << "  ";
        cout << endl;
    }       
}

