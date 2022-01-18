#include"rapidxml.hpp"
//#include"rapidxml_iterators.hpp"
#include"rapidxml_print.hpp"
#include"rapidxml_utils.hpp"
#include"iostream"
#include<fstream>
#include<string.h>
#include<iomanip>
#include<algorithm>
#include<math.h>
#include<queue>


//#define M_PI 3.14159265
using namespace rapidxml;
using namespace std;

struct edgeElement{
    int nodeIndex;
    string wayid;
    double distance;
};

struct node{
    long long int nodeid;
    string name;
    double latitude;
    double longitude;
};

struct heapElement{
    int index;
    double distance;
};

struct nodeDistance{
    struct node n;
    double distance;     
};
struct CompareDistances{
    bool operator()(heapElement const& n1, heapElement const& n2){
        return n1.distance>n2.distance;
    }
};
struct way{
    string wayid;
    //string name;
    vector<long long int> wayList;
};

struct parent{
    int index;
    string wayid;
};

//toRadians and distance is standard code implemented from online sources to calculate distance
double toRadians(const double degree)
{
    // cmath library in C++
    // defines the constant
    // M_PI as the value of
    // pi accurate to 1e-30
    double one_deg = (M_PI) / 180;
    return (one_deg * degree);
} 
double distance(struct node n1, struct node n2)
{
    // Convert the latitudes
    // and longitudes
    // from degree to radians.
    double lat1, long1, lat2, long2;
    
    lat1 = toRadians(n1.latitude);
    long1 = toRadians(n1.longitude);
    lat2 = toRadians(n2.latitude);
    long2 = toRadians(n2.longitude);
     
    // Haversine Formula
    double dlong = long2 - long1;
    double dlat = lat2 - lat1;
 
    double ans = pow(sin(dlat / 2), 2) +
                          cos(lat1) * cos(lat2) *
                          pow(sin(dlong / 2), 2);
 
    ans = 2 * asin(sqrt(ans));
 
    // Radius of Earth in
    // Kilometers, R = 6371
    // Use R = 3956 for miles
    double R = 6371;
     
    // Calculate the result
    ans = ans * R;
 
    return ans;
}

int NumberNodes = 0;
int NumberWays = 0;

bool compareNodes(struct node n1, struct node n2){
    return(n1.nodeid<n2.nodeid);
}

bool compareDistances(struct nodeDistance n1, struct nodeDistance n2){
    return (n1.distance<n2.distance);
}

void Count(){
    file<> xmlFile("map.osm");
    xml_document<>doc;
    doc.parse<0>(xmlFile.data());    
    xml_node<> *node = doc.first_node();    
    node = node->first_node();    
    node = node->next_sibling();        
    
    while(node){
        if (!strcmp(node->name(),"node")){
            NumberNodes++;                        
        }
        node = node->next_sibling();
    }   
    
    node = doc.first_node();
    node = node->first_node();
    
    while(node){
        if (!strcmp(node->name(),"way")){            
            NumberWays++;
        }
        node = node->next_sibling();
    }    
    
    
}

void Parsing(struct node nodes[], struct way ways[]){
    file<> xmlFile("map.osm");
    xml_document<>doc;
    doc.parse<0>(xmlFile.data());    
    xml_node<> *node = doc.first_node();    
    node = node->first_node();    
    node = node->next_sibling();  
   
    ofstream fout;
    fout.open("nodes.txt");
    
    int i = 0;
    while(node){
        if (!strcmp(node->name(),"node")){
            //NumberNodes++;
            //cout<<"node added";
            nodes[i].name = "Unknown";
            fout<<node->name()<<"\n";
            for (xml_attribute<> *attr = node->first_attribute(); attr; attr = attr->next_attribute())
            {
                fout<<" "<<attr->name()<<": ";
                fout<<" "<<attr->value();
                fout<<"\n";
                if (!strcmp(attr->name(),"id")){
                    nodes[i].nodeid = stol(attr->value());
                }
                if(!strcmp(attr->name(), "lat")){
                    nodes[i].latitude = stod(attr->value());
                }
                if(!strcmp(attr->name(), "lon")){
                    nodes[i].longitude = stod(attr->value());
                }

            }
            xml_node<> *child = node->first_node();
            fout<<"\ndisplaying child values\n";
            while(child){
                
                for (xml_attribute<> *attr = child->first_attribute(); attr; attr = attr->next_attribute())
                {
                    fout<<" "<<attr->name()<<": ";
                    fout<<" "<<attr->value();
                    if(!strcmp(attr->value(),"name")){
                        xml_attribute<> *attr1 = attr->next_attribute();
                        nodes[i].name = attr1->value();
                    }
                    
                }
                fout<<"\n";
                child = child->next_sibling();
            }
            i++;

        }
        node = node->next_sibling();
    }
    fout.close();
    //cout<<"Number of nodes is: "<<NumberNodes;
    fout.open("ways.txt");
    node = doc.first_node();
    node = node->first_node();
    cout<<node->name();
    i= 0;
    int NumberEdges = 0;
    while(node){
        if (!strcmp(node->name(),"way")){
            //cout<<"way added";
            //NumberWays++;
            fout<<node->name()<<"\n";
            fout<<"way";
            for (xml_attribute<> *attr = node->first_attribute(); attr; attr = attr->next_attribute())
            {
                fout<<" "<<attr->name()<<": ";
                fout<<" "<<attr->value();
                fout<<"\n";
                if(!strcmp(attr->name(), "id")){
                    ways[i].wayid = attr ->value();
                }
            }
            xml_node<> *nd = node->first_node();
            int j= 0;
            while(nd){

                if(!strcmp(nd->name(),"nd")){
                    j++;
                    xml_attribute<> *attr = nd->first_attribute();

                    fout<<" node id: "<<attr->value();
                    ways[i].wayList.push_back(stol(attr->value()));
                }
                nd = nd->next_sibling();
            }
            NumberEdges+=j;
            i++;

        }
        node = node->next_sibling();
    }
    fout.close();
    //cout<<"number of Edges discorvered is: "<<NumberEdges;
    

}

int findIndex(struct node nodes[], long long int nodeId){
    int first = 0; int last = NumberNodes-1;
    int mid;    
    while(first<=last){
        mid = first+(last-first)/2;
        if (nodes[mid].nodeid==nodeId){
            return mid;
        }        
        else if (nodes[mid].nodeid<nodeId){
            first = mid+1;
        }
        else{
            last = mid-1;
        }
    }
    return -1;
}

void PrintClosestNodes(struct node nodes[], long long int id, int k){
    int index;
    for(int i = 0; i<NumberNodes; i++){
        if(nodes[i].nodeid == id)
        index = i;
    }
    struct nodeDistance nodeDistances[NumberNodes];
    for(int i = 0; i<NumberNodes; i++){
        
        nodeDistances[i].n = nodes[i]; 
        nodeDistances[i].distance = distance(nodes[i], nodes[index]); 
    }
    ofstream fout;
    fout.open("distances.txt");
    for(int i = 0; i<NumberNodes; i++){
        fout<<nodeDistances[i].n.nodeid<<"--> distance: "<<nodeDistances[i].distance<<'\n';
    }
    fout.close();

    sort(nodeDistances, nodeDistances+NumberNodes, compareDistances);
    cout<<"\nthe k closest nodes are: \n";
    for(int i = 1; i<=k; i++){
        cout<<i<<"--> id: "<<nodeDistances[i].n.nodeid<<" distance: "<<nodeDistances[i].distance<<"\n";
    }
}

void search(struct node nodes[], string c){
    cout<<"\ndetails of any nodes found, if they exist: \n";   
    
    int j= 0;
    for(int i = 0; i<NumberNodes; i++){
        string name = nodes[i].name;
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        transform(c.begin(), c.end(), c.begin(), ::tolower);
        if(name.find(c)!= string::npos){
            j++;
            cout<<j<<"--> id: "<<nodes[i].nodeid<<" name: "<<nodes[i].name<<"\n";
        }
    }
}

void addEdge(vector <edgeElement> adj[], struct node nodes[], long long int nodeId1, long long int nodeId2, string wayid){
    int index1 = findIndex(nodes, nodeId1);
    int index2 = findIndex(nodes, nodeId2);
    edgeElement element1;
    element1.nodeIndex = index1;
    element1.distance = distance(nodes[index1], nodes[index2]);
    element1.wayid = wayid;
    edgeElement element2;
    element2.nodeIndex = index2;
    element2.distance = element1.distance;
    element2.wayid = wayid;
    adj[index1].push_back(element2);
    adj[index2].push_back(element1);

}

void makeGraph(vector <edgeElement> adj[], struct node nodes[], struct way ways[]){    
    for(int i = 0; i<NumberWays; i++){
        for(int j = 1; j<ways[i].wayList.size(); j++){

            addEdge(adj, nodes, ways[i].wayList[j], ways[i].wayList[j-1], ways[i].wayid);
        }
    }
}

void printGraph(vector <edgeElement> adj[]){
    ofstream fout;
    fout.open("graph.txt");
    for(int i = 0; i<NumberNodes; i++){
        fout<<i<<"connected to--> ";
        for(auto it: adj[i]){
            fout<<"     "<<it.nodeIndex<<" via: "<<it.wayid<<"distance: "<<it.distance;
        }
        fout<<"\n";
    }
    fout.close();
}

//minheap and dijkstra implementations were done with online references, with modifications
//made for my specific use case

struct MinHeapNode
{
    int  v;
    double dist;
};
 
// Structure to represent a min heap
struct MinHeap
{     
    // Number of heap nodes present currently
    int size;    
   
    // Capacity of min heap
    int capacity; 
   
    // This is needed for decreaseKey()
    int *pos;   
    struct MinHeapNode **array;
};
 
// A utility function to create a
// new Min Heap Node
struct MinHeapNode* newMinHeapNode(int v, int dist)
{
    struct MinHeapNode* minHeapNode =
           (struct MinHeapNode*)
      malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}
 
// A utility function to create a Min Heap
struct MinHeap* createMinHeap(int capacity)
{
    struct MinHeap* minHeap = (struct MinHeap*) malloc(sizeof(struct MinHeap));
    minHeap->pos = (int *)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array =(struct MinHeapNode**)malloc(capacity * sizeof(struct MinHeapNode*));
    return minHeap;
}
 
// A utility function to swap two
// nodes of min heap.
// Needed for min heapify
void swapMinHeapNode(struct MinHeapNode** a, struct MinHeapNode** b)
{
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}
 
// A standard function to
// heapify at given idx
// This function also updates
// position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(struct MinHeap* minHeap, int idx)
{
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;
 
    if (left < minHeap->size && minHeap->array[left]->dist < minHeap->array[smallest]->dist )
      smallest = left;
 
    if (right < minHeap->size && minHeap->array[right]->dist < minHeap->array[smallest]->dist )
      smallest = right;
 
    if (smallest != idx)
    {
        // The nodes to be swapped in min heap
        MinHeapNode *smallestNode = minHeap->array[smallest];
        MinHeapNode *idxNode = minHeap->array[idx];
 
        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;
 
        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);
 
        minHeapify(minHeap, smallest);
    }
}

int isEmpty(struct MinHeap* minHeap)
{
    return minHeap->size == 0;
}
 
// Standard function to extract
// minimum node from heap
struct MinHeapNode* extractMin(struct MinHeap* minHeap)
{
    if (isEmpty(minHeap))
        return NULL;
 
    // Store the root node
    struct MinHeapNode* root = minHeap->array[0];
 
    // Replace root node with last node
    struct MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;
 
    // Update position of last node
    minHeap->pos[root->v] = minHeap->size-1;
    minHeap->pos[lastNode->v] = 0;
 
    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);
 
    return root;
}
 
// Function to decreasy dist value
// of a given vertex v. This function
// uses pos[] of min heap to get the
// current index of node in min heap
void decreaseKey(struct MinHeap* minHeap, int v, int dist)
{
    // Get the index of v in  heap array
    int i = minHeap->pos[v];
 
    // Get the node and update its dist value
    minHeap->array[i]->dist = dist;
 
    // Travel up while the complete
    // tree is not hepified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist)
    {
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] = (i-1)/2;
        minHeap->pos[minHeap->array[(i-1)/2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);
 
        // move to parent index
        i = (i - 1) / 2;
    }
}
 
// A utility function to check if a given vertex
// 'v' is in min heap or not
bool isInMinHeap(struct MinHeap *minHeap, int v)
{
   if (minHeap->pos[v] < minHeap->size)
     return true;
   return false;
}
 
// A utility function used to print the solution
void printArr(double dist[], int n)
{
    ofstream fout;
    fout.open("dijkstra.txt");
    fout<<"Vertex   Distance from Source\n";
    for (int i = 0; i < n; ++i){
        //printf("%d \t\t %d\n", i, dist[i]);
        fout<<i<<"      "<<dist[i]<<"\n";
    }
    fout.close();
}
 
void printPath(struct node nodes[], struct parent parents[], int src, int dest){
    vector<parent> path;
    int j = 0;
    parent p;
    p.index = dest;
    p.wayid = "none";
    
    while(1){
        
        path.push_back(p);
        j++;
        if(p.index == src){
            break;
        }
        p = parents[p.index];
    }
    ofstream fout;
    fout.open("path.txt");
    cout<<"\ninstructions for shortest path can be found in path.txt: \n";
    for(int i = j-1; i>=1; i--){
        fout<<"go to: "<<nodes[path[i-1].index].nodeid<<" via wayid: "<<path[i].wayid<<"\n";
    }
    
}

// The main function that calculates
// distances of shortest paths from src to all
// vertices. It is a O(ELogV) function
void dijkstra(vector <edgeElement> adj[], int src, int dest, struct node nodes[])
{    
    int V = NumberNodes;  
    // dist values used to pick
    // minimum weight edge in cut
    double dist[V];  
    struct parent parents[V];  
 
    // minHeap represents set E
    struct MinHeap* minHeap = createMinHeap(V);
 
    // Initialize min heap with all
    // vertices. dist value of all vertices
    for (int v = 0; v < V; ++v)
    {
        dist[v] = 500000;
        minHeap->array[v] = newMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }
 
    // Make dist value of src vertex
    // as 0 so that it is extracted first
    minHeap->array[src] = newMinHeapNode(src, dist[src]);
    minHeap->pos[src]   = src;
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);
 
    // Initially size of min heap is equal to V
    minHeap->size = V;
 
    // In the followin loop,
    // min heap contains all nodes
    // whose shortest distance
    // is not yet finalized.
    while (!isEmpty(minHeap))
    {
        // Extract the vertex with
        // minimum distance value
        struct MinHeapNode* minHeapNode = extractMin(minHeap);
       
        // Store the extracted vertex number
        int u = minHeapNode->v;
 
        // Traverse through all adjacent
        // vertices of u (the extracted
        // vertex) and update
        // their distance values
        //struct AdjListNode* pCrawl =
        //             graph->array[u].head;
        for(auto it: adj[u])
        {
            int v = it.nodeIndex;
 
            // If shortest distance to v is
            // not finalized yet, and distance to v
            // through u is less than its
            // previously calculated distance
            if (isInMinHeap(minHeap, v) && dist[u] != 500000 && it.distance + dist[u] < dist[v])
            {
                dist[v] = dist[u] + it.distance;
                struct parent p;
                p.index = u;
                p.wayid = it.wayid;
                parents[v] = p;
 
                // update distance
                // value in min heap also
                decreaseKey(minHeap, v, dist[v]);
            }            
        }
    }

    cout<<"\nthe shortest distance between the 2 nodes is: "<<dist[dest];

    printPath(nodes, parents, src, dest);
 
    // print the calculated shortest distances
    printArr(dist, V);
}

int main(){
    Count();
    
    struct node nodes[NumberNodes];
    struct way ways[NumberWays];
    vector<edgeElement> adj[NumberNodes];
    Parsing(nodes, ways);
    // long long int i1 = findIndex(nodes, 3258023129);
    // long long int i2 = findIndex(nodes, 9236071057);

    // cout<<"\n\n distance: "<<distance(nodes[i1], nodes[i2])<<"\n\n";
    
    cout<<"\n\nenter 1 to find the total number of nodes and ways discovered in the given osm file\n";
    cout<<"enter 2 to search for a node by typing in a substring\n";
    cout<<"enter 3 to find the k closest nodes to a given node\n";
    cout<<"enter 4 to find the least distance between any 2 nodes\n";
    int a;
    cin>>a;
    int index1, index2;
    
    string s;
    switch(a){
        case 1:
            cout<<"\nThe number of nodes discovered is: "<<NumberNodes<<"\n";
            cout<<"The number of ways discovered is: "<<NumberWays<<"\n";
            break;
        case 2:   
            cout<<"\nenter the input string here: ";         
            cin>>s;
            search(nodes, s);
            break;
        case 3:
            cout<<"\n"<<"enter the nodeid you wish to view around: ";
            long long int id;
            cin>>id;
            cout<<"\n\nenter the number of closest nodes you wish to display (k): ";
            int k;
            cin>>k;
            PrintClosestNodes(nodes, id, k);
            break;
        case 4:
            makeGraph(adj, nodes, ways);
            printGraph(adj);
            
            cout<<"enter the 2 node ids between which you wish to view the shortest distance: \n";
            long long int id1;
            long long int id2;
            cin>>id1>>id2;
            index1 = findIndex(nodes, id1);
            index2 = findIndex(nodes, id2);
            dijkstra(adj, index1, index2, nodes);
            break;
        default: 
            cout<<"wrong input";
    }
}
