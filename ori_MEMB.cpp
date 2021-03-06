////////read .txt; omit the first 6 lines; from 0 to N-1;///////////////

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>

using namespace std;

const int MAX=12501;
const int MAXR=10;
const int MINR=1;
int N=0;
int L=0;//num of links
float s,e;
struct network
{
    int b;//unburnt--->0; burnt--->1
    vector<int> neighbor;
    int geodis[MAX];
    int excmas; // calculate "excluded mass" of each node --- once
    int cendis; // "central distance" 
    int boxid;  // uncovered=-1;covered=0; centernode=1,2,3...
};
network node[MAX];

int readnet(network *);//return num of links
void printnet(network *, int);

int excmas(network *, int); // return node with max excluded mass
bool complete(network *); // determine whether all nodes are covered
int MEMB(network *, int); // return Nb
float boxcr(network *,int);//return lB* --- averaged real lB
float module(network *,int);//calculate Modularity
float moduleN(network *,int);//calculate Modularity_Newman's defination
void khub(network *, int);
float aveout(network * node, int Nb);

void burning(network *);//calculate node[i].geodis[j]; return averaged path

int main()
{
    srand(time(NULL));

    ofstream outFile;
    outFile.open("MEMB.txt");

    outFile << "#lb" << "\t" << "lb*" << "\t" << "M" << "\t" << "M_Newman" << "\t" << "Nb" << "\t" <<  "s" << "\t" << "e" << "\t" << "Lout" <<endl;

    //ofstream pajek;
    //char Filename[40]; 

    L = readnet(node); // N=# of nodes
    
    burning(node);
    cout << "burning Done!" << endl;

    int Nb;
    float lb;
    float Mb;
    float Mb_Newman;
    float Lout;
    for (int Rb=MINR; Rb<=MAXR; Rb++)
    {
        Nb=MEMB(node,Rb);
        lb=boxcr(node,Nb);// lB*
        Mb=module(node,Nb);// Modularity
        Mb_Newman=moduleN(node,Nb);// Modularity
        //outFile << 2*Rb+1 << "\t" << lb << "\t" << Mb << "\t" << Nb << endl;
        khub(node,Nb);
        Lout=aveout(node,Nb);
        //cout << "s= " << s << endl;
        //cout << "e= " << e << endl;
        outFile << 2*Rb+1 << "\t" << lb << "\t"<< Mb <<"\t" << Mb_Newman << "\t" << Nb << "\t" << s << "\t" << e << "\t" << Lout << endl;
    }
    outFile.close();
    return 0;
}

void khub(network *, int Nb)
{
    ofstream koutput;
    koutput.open("khub.txt");
    koutput << "#khub" << "\t" << "k" << "\t" << "hub_id" << "\t" << "nh"<< endl;

    int k[Nb];
    int khub[Nb];
    int hub_id[Nb];
    int n_hub[Nb];
    for (int i=0;i<Nb;i++)
    {
        k[i]=0;
        khub[i]=0;
        n_hub[i]=0;
    }
    for (int i=0;i<N;i++)
        for (int j=i+1;j<N;j++)
            if ((node[i].geodis[j]==1)&&(node[i].boxid!=node[j].boxid))
                {
                     k[node[i].boxid-1]++;
                     k[node[j].boxid-1]++;
                }
    for (int i=0;i<N;i++)
        if (khub[node[i].boxid-1]<int(node[i].neighbor.size()))
        {
            khub[node[i].boxid-1]=node[i].neighbor.size();
            hub_id[node[i].boxid-1]=i;
        }
    for (int i=0; i<Nb; i++)
        for (int j=i+1; j<Nb; j++)
            if ((node[hub_id[i]].geodis[hub_id[j]]==1))
            {
                n_hub[i]++;
                n_hub[j]++;
            }       

    for (int i=0;i<Nb;i++)
        koutput << khub[i] << "\t" << k[i] << "\t" << hub_id[i] << "\t" << n_hub[i] << endl;

    koutput.close(); 

    s=0.0;
    e=0.0;
    for (int i=0;i<Nb;i++){
        s+= k[i]/float(khub[i]);
        e+= n_hub[i]/float(k[i]);
    }
    s/=Nb;
    e/=Nb;
    //cout << "s= " << s << endl;
    //cout << "e= " << e << endl;

}

//---calculate Modularity

float moduleN(network * node, int Nb)
{
    int sumd[Nb];
    float in[Nb];
    for (int i=0;i<Nb;i++)
    {
        sumd[i]=0;
        in[i]=0;
    }
    for (int i=0;i<N;i++)
        for (int j=i+1;j<N;j++)
            if (node[i].geodis[j]==1)
                if (node[i].boxid==node[j].boxid)
                    in[node[i].boxid-1]++;
    for (int i=0; i<N; i++)
        sumd[node[i].boxid-1] += node[i].neighbor.size();

    float m=0;
    for (int i=0;i<Nb;i++)
        m += in[i]/L-sumd[i]*sumd[i]/4/L/L;
    m=m/Nb;
    cout << m << endl;
    return m;
}

// burning ---> calculate node[i].geodis[j]

void burning(network * node)
{
    int i,t,c,g,h;
    //int ldis=0;//-------Longest geodistance
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
        g=0;
        while (active.size()) 
        {
            g++;
            if (g>2*MAXR) break;
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
        }
        //g--;
        //if (g>ldis) ldis = g;

    }

    /*
    int sum_dis=0;
    for (int i=0;i<N;i++)
        for (int j=i+1;j<N;j++)
            sum_dis+= node[i].geodis[j];
    float L=float(sum_dis)/(N*(N-1)/2);
    return L;
    */
}


//---calculate Modularity

float aveout(network * node, int Nb)
{
    float out[Nb];
    for (int i=0;i<Nb;i++)
    {
        out[i]=0;
    }
    for (int i=0;i<N;i++)
        for (int j=i+1;j<N;j++)
            if (node[i].geodis[j]==1)
                if (node[i].boxid!=node[j].boxid)
                {
                     out[node[i].boxid-1]++;
                     out[node[j].boxid-1]++;
                }
    float ave_out=0.0;
    for (int i=0;i<Nb;i++)
        ave_out+=out[i];
    ave_out/=Nb;
    return ave_out;
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
    /*float ave_out=0.0;
    for (int i=0;i<Nb;i++)
        ave_out+=out[i];
    ave_out/=Nb;
    cout << "ave_out= " << ave_out << endl;
    */
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

int MEMB(network * node, int Rb)
{
    for (int i=0;i<N;i++)
    {
        node[i].boxid = -1; // mark all nodes as uncovered and noncenters
        node[i].excmas = 0;
        node[i].cendis = 0;
    }

    //---calculation of the box centers
    int Nb=1;
    int p;
    while (!complete(node))
    {
        p=excmas(node,Rb);
        //cout << "center_id= " << p << endl;
        node[p].boxid=Nb;
        for (int j=0;j<N;j++)
            if ((node[p].geodis[j]<=Rb)&&(node[p].geodis[j]>0)&&(node[j].boxid<0))
                node[j].boxid=0;
        Nb++;
        //cout << "." << endl;
    }
    Nb=Nb-1;
    //cout << "Nb= " << Nb << endl;
    
    ofstream boxcenter;
    boxcenter.open("boxcenter.txt");
    
    for (int i=0;i<N;i++)
        if (node[i].boxid>0)
            boxcenter << i << "\t" << node[i].boxid << endl;

    boxcenter.close();

    
    //---calculation of the box centers---END

    
    //---assign each node to its nearest center
    int mini;
    for (int i=0;i<N;i++)
        if (node[i].boxid==0)
        {
            mini=Rb;
            for (int j=0;j<N;j++)
                if ((node[i].geodis[j]<mini)&&(node[i].geodis[j]>0)&&(node[j].boxid>0))
                    mini=node[i].geodis[j];
            node[i].cendis=mini;
        }
    /*
    for (int i=3;i<=3;i++)
        cout << i << "'s cendis is " << node[i].cendis << endl;
    */
    for (int t=1;t<=Rb;t++)
        for (int i=0;i<N;i++)
            if ((node[i].boxid==0)&&(node[i].cendis==t))
            {
                int s=node[i].neighbor.size();
                int u[s];
                int e=0;
                for (int j=0;j<s;j++)
                {
                    int w=node[i].neighbor[j];
                    if ((node[w].cendis<node[i].cendis)&&(node[w].boxid>0))
                    {
                        u[e]=node[w].boxid;
                        e++;
                    }
                }
                if (e==0)
                    cout << "e=0"<< endl;
                else
                {
                    int z=rand()%e;
                    node[i].boxid=u[z];
                }
                //if (i==3) cout << "e= " << e << endl;

            }
            
    ofstream outFile;
    outFile.open("boxid.txt");
    
    for (int i=0;i<N;i++)
        outFile << i << "\t" << node[i].boxid << endl;

    outFile.close();
    
    //---assign each node to its nearest center---END

    //printnet(newnode, Nb);
    
    return Nb;

}

// --- determine whether all nodes are covered

bool complete(network * node)
{
    bool flag=true;
    for (int i=0;i<N;i++)
        if (node[i].boxid<0)
            flag=false;
    //cout << "complete check done! " << endl;
    return flag;
}

// --- calculate "excluded mass"

int excmas(network * node, int Rb)
{
    for (int i=0;i<N;i++){
        if (node[i].boxid==-1)
            node[i].excmas = 1;
        else node[i].excmas = 0;
    }

    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
            if ((node[i].geodis[j]<=Rb)&&(node[i].geodis[j]>0)&&(node[j].boxid==-1))
                node[i].excmas++;
    /*
    for (int i=0;i<N;i++)
        cout << i << "'s excmas = " << node[i].excmas << endl;
    */
    int mem=0;
    int mem_id=0;
    for (int i=0;i<N;i++)
        if ((node[i].excmas>mem)&&(node[i].boxid<=0))
        {
            mem=node[i].excmas;
            mem_id=i;
        }
    //cout << mem_id << endl;

    return mem_id;
}

// --- read network & get N

int readnet(network * node)    
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
        node[i].boxid = -1; // mark all nodes as uncovered and noncenters
        node[i].excmas = 0;
        node[i].cendis = 0;
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
    inFile.close();

    for (int i=0;i<MAX;i++)
        if (node[i].neighbor.size())
            N++;
    cout <<  "# of nodes: " << N << endl;

    int sum_k=0;
    for (int i=0;i<N;i++)
        sum_k+=node[i].neighbor.size();
    sum_k/=2;
    cout << "# of links : " << sum_k << endl;

    return sum_k;

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

