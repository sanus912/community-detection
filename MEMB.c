#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct node
{
        int data;
        struct node *next;
} *NODEPTR;

int *ivector(long nl, long nh);
void add_node(int i, struct node** headRef);
void clear_lists(int cluster_size);
void free_ivector(int *v, long nl, long nh);

struct node **neigh_list;

int *neighbors;

struct node **re_neigh_list;


int main(argc,argv)
int argc;
char *argv[];
{
	int i,j,cluster_size,boxno,l_b,paint_box();
	char name[50];
	int are_neighbors(),*color;
	int unboxed,new_size,old_size,id1,id2,stage,*re_neighbors,*dist_to_hub,d;
	struct node *current,*cur2,*hubs_list;
	void locate_hubs();
	FILE *fp;
	
	if(argc==4)
	{
		cluster_size=atoi(argv[2]);
		l_b=atoi(argv[3]);
	}
	else
        {
                printf("Usage: %s filename number_of_nodes r_B\n",argv[0]);
                return(-1);
        }

	/* Initialize random number generator */
	srand48(0);

	/* Read the network from file */
	/* Network has to be two columns with the links. Nodes should be numbered from 0 to N-1 */
	/* All node indexes should be present */
	/* Each link should appear once (e.g. "A B" and "B A" are not allowed */
        fp=fopen(argv[1],"r");
        neighbors=ivector(0,cluster_size-1);
        for(i=0;i<cluster_size;i++) neighbors[i]=0;
        /* Allocate neighbor lists memory for each node */
        neigh_list=(struct node **)malloc(cluster_size*sizeof(struct node));
        for(i=0;i<cluster_size;i++)
        {
                neigh_list[i]=(struct node *)malloc(sizeof(struct node));
                neigh_list[i]=NULL;
        }

        while(!feof(fp))
        {
                fscanf(fp,"%d\t%d\n",&i,&j);
                neighbors[i]++;
                neighbors[j]++;
		add_node(i,&neigh_list[j]);
		add_node(j,&neigh_list[i]);
        }
        fclose(fp);
	
	hubs_list=(struct node *)malloc(sizeof(struct node));
	hubs_list=NULL;

	color=ivector(0,cluster_size-1);
	dist_to_hub=ivector(0,cluster_size-1);
	
	old_size=cluster_size;
	new_size=cluster_size;
	stage=1;
        //while(new_size>1)
	//{
		for(i=0;i<old_size;i++)
		{
			color[i]=0;
			dist_to_hub[i]=-1;
		}
		unboxed=old_size;
		boxno=1;
		locate_hubs(old_size,l_b,&hubs_list);
		
		current=hubs_list;
		while(current!=NULL)
		{
			color[current->data]=boxno;
			dist_to_hub[current->data]=0;
			unboxed--;
			current=current->next;
			boxno++;
		}
		current=hubs_list;
		while(current!=NULL)
		{
			cur2=neigh_list[current->data];
			while(cur2!=NULL)
			{
				if(color[cur2->data]==0)
				{
					color[cur2->data]=color[current->data];
					dist_to_hub[cur2->data]=1;
					unboxed--;
				}
				else if(dist_to_hub[cur2->data]==1 && drand48()<0.5) color[cur2->data]=color[current->data];
				cur2=cur2->next;
			}
			current=current->next;
		}
		i=0;
		d=2;
		while(unboxed>0)
		{
			if(color[i]==0)
			{
				current=neigh_list[i];
				while(current!=NULL && color[i]==0)
				{
					if(dist_to_hub[current->data]==d-1)
					{
						color[i]=color[current->data];
						dist_to_hub[i]=d;
						unboxed--;
					}
					else current=current->next;
				}
			}
			i=(i+1)%old_size;
			if(i==0) d++;
		}

		new_size=boxno-1;
		if(stage==1)
		{
			re_neighbors=ivector(0,new_size-1);
			/* Allocate neighbor lists memory for each node */
        		re_neigh_list=(struct node **)malloc(new_size*sizeof(struct node));
        		for(i=0;i<new_size;i++)
        		{
	        		re_neigh_list[i]=(struct node *)malloc(sizeof(struct node));
        			re_neigh_list[i]=NULL;
        		}
        	}
        	for(i=0;i<new_size;i++) re_neighbors[i]=0;
		for(i=0;i<old_size;i++)
		{
			id1=color[i]-1;
			current=neigh_list[i];
			while(current!=NULL)
			{
				id2=color[current->data]-1;
				if(id1!=id2)
				{
					if(!are_neighbors(id1,id2))
					{
						add_node(id1,&re_neigh_list[id2]);
						add_node(id2,&re_neigh_list[id1]);
						re_neighbors[id1]++;
						re_neighbors[id2]++;
					}
				}
				current=current->next;
			}
		}
		for(i=0;i<old_size;i++)
        	{
        		current=neigh_list[i];
        		while(current!=NULL)
        		{
                		cur2=current->next;
                		free(current);
                		current=cur2;
			}
			neigh_list[i]=NULL;
		}
		if(new_size>1)
		{
			sprintf(name,"k_vs_khub_stage%d_rb_%d.dat",stage-1,l_b);
			fp=fopen(name,"w");
			current=hubs_list;
			while(current!=NULL)
			{
				fprintf(fp,"%d\t%d\t%d\n",neighbors[current->data],re_neighbors[color[current->data]-1],current->data);
				current=current->next;
			}
			fclose(fp);
		}
		if(new_size>1)
		{
			sprintf(name,"boxes_stage%d_rb_%d.dat",stage-1,l_b);
			fp=fopen(name,"w");
			for(i=0;i<old_size;i++) fprintf(fp,"%d\t%d\n",i,color[i]);
			fclose(fp);
		}
		for(i=0;i<new_size;i++)
		{
			neighbors[i]=re_neighbors[i];
			neigh_list[i]=re_neigh_list[i];
			current=re_neigh_list[i];
			re_neigh_list[i]=NULL;
		}
		old_size=new_size;
		
		if(new_size>1)
		{
			sprintf(name,"net_stage%d_rb_%d.dat",stage,l_b);
			fp=fopen(name,"w");
			for(i=0;i<new_size;i++)
			{
				current=neigh_list[i];
				while(current!=NULL)
				{
					if(i<j) fprintf(fp,"%d\t%d\n",i,current->data);
					current=current->next;
				}
			}
		}
		
        	stage++;
        	current=hubs_list;
		while(current!=NULL)
		{
			cur2=current->next;
			free(current);
			current=cur2;
		}
		hubs_list=NULL;
	//}
	free_ivector(re_neighbors,0,cluster_size-1);
	free_ivector(color,0,cluster_size-1);
	free_ivector(dist_to_hub,0,cluster_size-1);
	
	clear_lists(cluster_size);
	
	return(0);
}

int paint_box(origin,l_b,size,boxno,color)
int origin,size,l_b,boxno,*color;
{
        int i,distance,n_inf,n,c2,total,found_new;
        int *just_inf,*infected,*visited;
        struct node *temp;

	visited=ivector(0,size-1);
        infected=ivector(0,size-1);
        just_inf=ivector(0,size-1);

	for(i=0;i<size;i++) visited[i]=0;
        n_inf=1;
        infected[0]=origin;
        visited[origin]=1;
        distance=1;
        found_new=0;
        do
        {
                total=0;
                for(i=0;i<n_inf;i++)
                {
                        n=infected[i];
                        temp=neigh_list[n];
                        while(temp!=NULL)
                        {
                                c2=temp->data;
                                if(visited[c2]==0)
                                {
                                        visited[c2]=1;
                                        just_inf[total]=c2;
                                        total++;
                                        if(color[c2]==0)
                                        {
                                        	color[c2]=boxno;
                                        	found_new++;
                                        }
                                }
                                temp=temp->next;
                        }
                }
                n_inf=total;

                for(i=0;i<total;i++) infected[i]=just_inf[i];
                distance++;
        } while(distance<=l_b);

	free_ivector(visited,0,size-1);
        free_ivector(just_inf,0,size-1);
        free_ivector(infected,0,size-1);

        return(found_new);
}

int are_neighbors(cand,index)
int cand,index;
{

        struct node *current;

        if(cand==index) return 1;

        current=re_neigh_list[index];
        while(current!=NULL)
        {
                if(cand==current->data) return 1;
                else current=current->next;
        }
        return 0;
}

void locate_hubs(size,l_b,hubs_list)
int size,l_b;
struct node *hubs_list;
{
        int i,distance,n_inf,n,c2,total,found_new,origin,hub,m_max,unassigned;
        int *just_inf,*infected,*visited,*mass,*assigned,*center_candidate,paint_box(),nhubs;
        struct node *temp;
        void add_bottom_node();

	mass=ivector(0,size-1);
	visited=ivector(0,size-1);
        infected=ivector(0,size-1);
        just_inf=ivector(0,size-1);
        assigned=ivector(0,size-1);
        center_candidate=ivector(0,size-1);

	for(i=0;i<size;i++)
	{
		assigned[i]=0;
		center_candidate[i]=1;
	}
	unassigned=size;
	nhubs=0;
	while(unassigned>0)
	{
		for(i=0;i<size;i++)
		{
			if(assigned[i]==0) mass[i]=1;
			else mass[i]=0;
		}
        	for(origin=0;origin<size;origin++)
        	{
			if(center_candidate[origin]==1)
			{
				for(i=0;i<size;i++) visited[i]=0;
        			n_inf=1;
        			infected[0]=origin;
        			visited[origin]=1;
        			distance=1;
        			found_new=0;
				do
				{
					total=0;
					for(i=0;i<n_inf;i++)
					{
						n=infected[i];
						temp=neigh_list[n];
						while(temp!=NULL)
						{
							c2=temp->data;
							if(visited[c2]==0)
							{
								visited[c2]=1;
								just_inf[total]=c2;
								total++;
								if(assigned[c2]==0) mass[origin]++;
                                			}
                                			temp=temp->next;
                        			}
                			}
                			n_inf=total;

                			for(i=0;i<total;i++) infected[i]=just_inf[i];
                			distance++;
        			} while(distance<=l_b);
				if(mass[origin]==0) center_candidate[origin]=0;
        		}
		}
	
		m_max=0;
		hub=-1;
		for(i=0;i<size;i++)
		{
			if(m_max<mass[i])
			{
				m_max=mass[i];
				hub=i;
			}
		}
		nhubs++;
		center_candidate[hub]=0;
		add_bottom_node(hub,hubs_list);

		if(assigned[hub]==0) assigned[hub]=hub+1;
		paint_box(hub,l_b,size,hub+1,assigned);
		unassigned-=m_max;
	}

	free_ivector(center_candidate,0,size-1);
	free_ivector(mass,0,size-1);
	free_ivector(visited,0,size-1);
        free_ivector(just_inf,0,size-1);
        free_ivector(infected,0,size-1);
        free_ivector(assigned,0,size-1);

	return;
}

void add_bottom_node(i,headRef)
int i;
struct node** headRef;
{
        struct node *newnode,*current;

        newnode = (struct node *)malloc(sizeof(struct node));
        newnode->data=i;
        newnode->next=NULL;
        if(*headRef==NULL) *headRef=newnode;
        else
        {
		current=*headRef;
		while(current->next!=NULL) current=current->next;
		current->next=newnode;
	}

        return;
}

void add_node(i,headRef)
int i;
struct node** headRef;
{ 
        struct node *newnode;
 
        newnode = (struct node *)malloc(sizeof(struct node));
        newnode->data=i;
        newnode->next=*headRef;
 
 	*headRef=newnode;
	
        return;
}

void clear_lists(cluster_size)
int cluster_size;
{
	struct node *current,*cur2;
	int i;
	
	for(i=0;i<cluster_size;i++)
	{
		current=neigh_list[i];
               	while(current!=NULL)
               	{
                       	cur2=current->next;
                       	free(current);
                       	current=cur2;
	       	}
       	       	neigh_list[i]=NULL;
	}
	free(neigh_list);
	free_ivector(neighbors,0,cluster_size-1);
	return;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;
        v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
        if (!v) {printf("allocation failure in ivector()\n");exit(-1);}
        return v-nl+1;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((char*) (v+nl-1));
}
