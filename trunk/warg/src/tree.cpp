#include "tree.h"
#include "recedge.h"
#include "slotallocator.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stack>
#include <limits>

using namespace std;
namespace weakarg {

    Tree::Tree(int n) {
      this->n=n;
      nodes=vector<Node*>(n+n-1);
      Node ** toCoal=(Node **)calloc(n,sizeof(Node*));
      int k=n;
      double time=0.0;
      for (int i=0;i<n;i++) {
          nodes[i]=SlotAllocator<Node>::GetSlotAllocator().Allocate();
          new (nodes[i]) Node();
          nodes[i]->setId(i);
          toCoal[i]=nodes[i];
        }
      while (k>1) {
          time-=2.0/(k*(k-1.0))*log(gsl_rng_uniform(rng));
          int a=(int)floor(gsl_rng_uniform(rng)*k);
          int b=a;
          while (b==a)
            b=(int)floor(gsl_rng_uniform(rng)*k);
          nodes[n+n-k]=SlotAllocator<Node>::GetSlotAllocator().Allocate();
          new (nodes[n+n-k]) Node();
          nodes[n+n-k]->setId(n+n-k);
          nodes[n+n-k]->setLeft(toCoal[a]);
          nodes[n+n-k]->setRight(toCoal[b]);
          nodes[n+n-k]->setAge(time);
          toCoal[a]->setFather(nodes[n+n-k]);
          toCoal[b]->setFather(nodes[n+n-k]);
          toCoal[a]=nodes[n+n-k];
          toCoal[b]=toCoal[k-1];
          k--;
        }
      root=nodes[n+n-2];
      free(toCoal);
      computeTTotal();
    }

    void Tree::makeFromNewick(std::string newick,bool forceages) {
      //Remove final ';'
      if (newick.at(newick.length()-1)==';')
        newick=newick.substr(0,newick.length()-1);
      vector<Node*> all;
      //Create tree
      root=SlotAllocator<Node>::GetSlotAllocator().Allocate();
      new (root) Node(&newick,NULL,&all);	// use in-place constructor
      n=(all.size()/2)+1;
      vector<double> dists;
      for (unsigned int i=0;i<all.size();i++) {
          if (all[i]->getLeft()==NULL)
            dists.push_back(0.0);
          else
            dists.push_back(all[i]->getLeft()->getAge());
        }
      for (unsigned int i=0;i<all.size();i++) {
          if (all[i]->getLeft()==NULL)
            all[i]->setAge(0.0);
          else
            all[i]->setAge(all[i]->getLeft()->getAge()+dists[i]);
        }
      vector<double> ages;
      for (unsigned int i=0;i<all.size();i++)
        ages.push_back(all[i]->getAge());
      //Create list of nodes ordered by age
      nodes=vector<Node*>(n+n-1);
      // This needs to be changed in order to accept newick files with labels not from 0:(n-1)
	// check whether the smallest id is 0 or some other number.
	int idbase = (std::numeric_limits<int>::max)();
	for (unsigned int i=0;i<all.size();i++) {
		if(all[i]->getId()<idbase && all[i]->getId() >= 0)
			idbase = all[i]->getId();
	}
	if(idbase > (int)all.size()||idbase<0)	idbase = 0;
	for (unsigned int i=0;i<all.size();i++) {
		if (all.at(i)->getId()!=-1){
			all.at(i)->setId(all.at(i)->getId()-idbase);
		}
	}
      for (unsigned int i=0;i<all.size();i++) {
          if (all.at(i)->getId()==-1) {
              int idnew=0;
              for (unsigned int j=0;j<all.size();j++)
                if (ages[j]<ages[i] || (ages[j]==ages[i] && j<i))
                  idnew++;
              all.at(i)->setId(idnew);
            }
          nodes.at(all.at(i)->getId())=all.at(i);
        }
      int incomcounts=0;
      for (unsigned int i=n;i<nodes.size();i++) {
          if (forceages && nodes.at(i)->getAge()<=nodes.at(i-1)->getAge()) {
              incomcounts++;
              nodes[i]->setAge(nodes[i-1]->getAge()+0.001);
            }
        }
      if (forceages && incomcounts>0) cerr<<"WARNING: Tree node ages are not strictly increasing.  Tree has been deterministically perturbed to ensure consistency with the coalescent model! ("<<incomcounts<<" Nodes affected)"<<endl;

      computeTTotal();
    }

    Tree::Tree(string newickfile,bool isFilename,bool forceages) {
      string newick;
      if (isFilename) {
          //Read Newick string from file
          ifstream newickFile;
          newickFile.open(newickfile.data());
          if (!newickFile) {
              cerr << "Can't open input file " << newickfile.data() << endl;
              exit(1);
            }
	  getline(newickFile,newick);
          newickFile.close();
          makeFromNewick(newick,forceages);
        } else makeFromNewick(newickfile,forceages);
	//cout<<"Using tree "<<root->newick()<<endl;
    }

    Tree::Tree(const Tree& t) {
    	root=NULL;
    	assign(t);
    }

    Tree::Tree(Tree *intree) {
      makeFromNewick(intree->newick(),false);
      // correct for any rounding
      for (int i=0;i<n;i++) {
          if (nodes.at(i)->getAge()!=intree->getNode(i)->getAge()) cout<<"node "<<i<<" dist="<<nodes.at(i)->getAge()-intree->getNode(i)->getAge()<<endl;
          nodes.at(i)->setAge(intree->getNode(i)->getAge());
        }
      ttotal=intree->getTTotal();
    }

    string Tree::newick(int p) const {
        return root->newick(p)+";";
      }

    string Tree::newickNoInternalLabels(int p) const {
        return root->newickNoInternalLabels(p)+";";
      }

    Tree::~Tree() {
      if (root)
        SlotAllocator<Node>::GetSlotAllocator().Free(root);
      root=NULL;
    }

    void Tree::assign(const Tree& t){
    	n = t.n;
		ttotal=t.ttotal;
		if(root!=NULL) SlotAllocator<Node>::GetSlotAllocator().Free(root);
    	root = t.root->cloneSubtree();
    	// put nodes in the node array at a position given by ID
    	nodes.resize(t.nodes.size());
    	stack<Node*> s;
    	s.push(root);
    	while(s.size()>0)
    	{
    		Node* cur = s.top();
    		s.pop();
    		// set values for the new node
    		nodes[cur->getId()] = cur;
    		if(cur->getLeft() != NULL)
    			s.push(cur->getLeft());
    		if(cur->getRight() != NULL)
    			s.push(cur->getRight());
    	}
    }

    double Tree::prior() const {
        double ret=0.0;
        for (int i=2;i<=n;i++)
          ret-= 0.5*i*(i-1)*(nodes[2*n-i]->getAge()-nodes[2*n-i-1]->getAge());
        return ret;
      }

    void Tree::computeTTotal() {
      ttotal=0.0;
      for (unsigned int i=n;i<nodes.size();i++)
        ttotal+=(nodes[i]->getAge()-nodes[i-1]->getAge())*(2*n-i);
    }

    int Tree::getPoint(double * dist, std::vector<int> * samplespace) const {
        int ans=-1;
        while (ans<0) {
            *dist=gsl_rng_uniform(rng)*ttotal;
            int i=0;
            while (1) {
                if (*dist>nodes[i]->getDist()) *dist-=nodes[i++]->getDist();
                else {
                    if (samplespace==NULL) ans = i;
                    else {
                        for (unsigned int j=0;j<samplespace->size();j++) if (samplespace->at(j)==i) {ans=i;j=samplespace->size();}
                      }
                    break;
                  }
              }
          }
        return(ans);
      }

    void Tree::orderNodes(double dist) {
      if (dist<0) for (int i=2*getN()-2;i>=getN();i++) orderNodes(i,dist);
      if (dist>0) for (int i=getN();i<2*getN()-1;i++) orderNodes(i,dist);
    }

    int Tree::orderNodes(int which, double dist) {
      int newloc=-1;
      // keep swapping which with its neighbour until the ordering is correct
      if (dist<0) {
          for (int i=which;i>=getN();i--) {
              if (getNode(i-1)->getAge()>getNode(i)->getAge())
                swapNode(i-1,i);
              else {
                  newloc=i;
                  break;
                }
            }
        } else {
          for (int i=which+1;i<2*getN()-1;i++) {
              if (getNode(i-1)->getAge()>getNode(i)->getAge())
                swapNode(i-1,i);
              else {
                  newloc=i-1;
                  break;
                }
            }
        }
      if (newloc==-1)
        newloc=2*getN()-2; // if this is the oldest node it should stay where it is
      return newloc;
    }

    void Tree::swapNode(int a, int b) {
		swap(nodes[a],nodes[b]);
      nodes[a]->setId(a);
      nodes[b]->setId(b);

    }

    void Tree::swapFather(int a,int b) {
      Node * fa=nodes[a]->getFather();
      Node * fb=nodes[b]->getFather();
      bool faleft=false,fbleft=false;
      if (fa!=NULL) {
          if (fa->getLeft()->getId()==a) {
              faleft=true;
            } else if (fa->getRight()->getId()==a) {
              faleft=false;
            } else {
              cout<<"Node fa "<<a<<" has father "<<fa->getId()<<endl;
              cout<<"But father has daughters "<<fa->getLeft()->getId()<<" and "<<fa->getRight()->getId()<<endl;
              throw("Daughter and father do not match!");
            }
        }
      if (fb!=NULL) {
          if (fb->getLeft()->getId()==b) {
              fbleft=true;
            } else if (fb->getRight()->getId()==b) {
              fbleft=false;
            } else {
              cout<<"Node fb "<<b<<" has father "<<fb->getId()<<endl;
              cout<<"But father has daughters "<<fb->getLeft()->getId()<<" and "<<fb->getRight()->getId()<<endl;
              throw("Daughter and father do not match!");
            }
        }
      if (faleft && fa!=NULL) fa->setLeft(nodes[b]);
      else if (fa!=NULL) fa->setRight(nodes[b]);
      if (fbleft && fb!=NULL) fb->setLeft(nodes[a]);
      else if (fb!=NULL) fb->setRight(nodes[a]);
      nodes[a]->setFather(fb);
      nodes[b]->setFather(fa);
      fa=nodes[a]->getFather();
      fb=nodes[b]->getFather();
      if (fa==NULL) root=nodes[a];
      if (fb==NULL) root=nodes[b];
    }

    void Tree::testNodeAges() const {
        for (int i=0;i<getN()*2-2;i++) {
            if (getNode(i)->getId()!=i) {
                cerr<<"Node "<<i<<" has ID "<<getNode(i)->getId()<<endl;
                throw "tree:ID labelling Broken";
              }
            if (getNode(i)->getAge()>getNode(i+1)->getAge()) {
                cerr<<"Age["<<i<<"]="<<getNode(i)->getAge()<<", Age["<<i+1<<"]="<<getNode(i+1)->getAge()<<endl;
                throw "tree:Age ordering Broken";
              }
          }
      }

    int Tree::getOldestReversedNode() {
      int onode=-1,n=getN();
      double curmaxage=0;
      for (int i=n;i<2*n-1;i++) {//all internal nodes...
          if (nodes[i]!=root) {//bar root node (at this stage they are not sorted)
              if (getNode(i)->getAge()>getNode(i)->getFather()->getAge()) {
                  if (getNode(i)->getAge()>curmaxage) {
                      curmaxage=getNode(i)->getAge();
                      onode=i;
                    }
                }
            }
        }
      return(onode);
    }

    std::vector<int> Tree::getAllChildren(int e) {
      vector<int> my;
      my.push_back(e);
      if (e<getN()) {return(my);} else {
          vector<int> rvec,lvec;
          lvec=getAllChildren(getNode(e)->getLeft()->getId());
          rvec=getAllChildren(getNode(e)->getRight()->getId());
          for (unsigned int i=0;i<rvec.size();i++) my.push_back(rvec[i]);
          for (unsigned int i=0;i<lvec.size();i++) my.push_back(lvec[i]);
        }
      return(my);
    }

    std::vector<int> Tree::getAllSampledSeqs(int e) {
      vector<int> my;
      if (e<getN()) {
	my.push_back(e);
 	return(my);
      }else {
          vector<int> rvec,lvec;
          lvec=getAllSampledSeqs(getNode(e)->getLeft()->getId());
          rvec=getAllSampledSeqs(getNode(e)->getRight()->getId());
          for (unsigned int i=0;i<rvec.size();i++) if (rvec[i]<getN()) my.push_back(rvec[i]);
          for (unsigned int i=0;i<lvec.size();i++) if (lvec[i]<getN()) my.push_back(lvec[i]);
        }
      return(my);
    }

    double Tree::tavare() const {
        double ret=log(0.5*(n-1));
        ret+=-0.5*getTTotal();
        ret+=(n-2)*log(1.0-exp(-0.5*getTTotal()));
        return ret;
      }

    vector<int> Tree::getMinEdgeList(vector<int> seqs) {
      vector<int> ret;
      while (seqs.size()>0) ret.push_back(getNextGroup(&seqs));
      return(ret);
    }


    int Tree::getNextGroup(vector<int> *seqs) {
      int par=seqs->at(0),thischild=seqs->at(0),ochild,retval=thischild;
      if ((int)seqs->size()==getN()) {
        seqs->clear();
	return(root->getId());
      }
      seqs->erase(seqs->begin());
      vector <int> tmplist;
      while (seqs->size()>0) {
          if (getNode(thischild)->getFather()==NULL) throw("Invalid node!");
          par=getNode(thischild)->getFather()->getId();
          ochild=otherChild(par,thischild);
          vector<int> listfound;
          if (isParentToOnly(ochild,*seqs,&listfound)) {
              for (unsigned int i=0;i<listfound.size();i++) for (unsigned int j=0;j<seqs->size();j++)  if (seqs->at(j)==listfound.at(i)) {
                      if (getNode(thischild)->getFather()->getDist()<0.00001) tmplist.push_back(seqs->at(j));
                      else retval=thischild;
                      seqs->erase(seqs->begin() + j);

                    }
            } else {// not a parent to only the seqs in our list
              for (unsigned int i=0;i<tmplist.size();i++) { seqs->push_back(tmplist[i]);}
              return(retval);
            }
          thischild=par;
        }
      if (tmplist.size()==0) retval=par;
    for (unsigned int i=0;i<tmplist.size();i++) { seqs->push_back(tmplist[i]);}
      return(retval);
    }

    bool Tree::isParentToOnly(int e, vector<int> seqs,vector<int> *which) {
      vector<int> ret;
      vector<int> allchildren=getAllSampledSeqs(e);
      for (unsigned int i=0;i<allchildren.size();i++) {
          bool found=false;
          for (unsigned int j=0;j<seqs.size();j++) {
              if (allchildren.at(i)==seqs.at(j)) {
                  found=true;
                  ret.push_back(seqs.at(j));
                  continue;
                }
            }
          if (!found) return(false);
        }
      if (which!=NULL) {
          which->clear();
          for (unsigned int i=0;i<ret.size();i++) which->push_back(ret[i]);
        }
      return(true);
    }

    Tree::Tree(Data * data) {
      this->n=data->getN();
      nodes=vector<Node*>(n+n-1);
      //First build distance matrix
      vector< vector<double> > dist(2*n,vector<double>(2*n,0.0));
      for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
            double d=0.0;
            double div=0.0;
            for (int k=0;k<data->getL();k++) {
                if (data->get(i,k)=='N' || data->get(j,k)=='N') continue;
                div++;
                if (data->get(i,k)!=data->get(j,k)) d++;
              }
            if(div==0){
            	cerr << "Warning! sequences " << i << "," << j << " have no aligned sites!\n";
            	div = 1.0;	// fudge it.
            }
            dist[i][j]=d/div+gsl_rng_uniform(rng)*1.0e-6;
            dist[j][i]=dist[i][j];
          }

      //Prepare to cluster
      vector<int> toCluster;
      vector<int> toClusterSize;
      for (int i=0;i<n;i++) {
          nodes[i]=SlotAllocator<Node>::GetSlotAllocator().Allocate();
          new (nodes[i]) Node();
          nodes[i]->setId(i);
          toCluster.push_back(i);
          toClusterSize.push_back(1);
        }
      int k=n;

      //Clustering procedure
      while (toCluster.size()>1) {
          //Find smallest distance
          double mini=0.0;
          int a=-1;
          int b=-1;
          for (unsigned int i=0;i<toCluster.size();i++) for (unsigned int j=i+1;j<toCluster.size();j++)
                if (dist[toCluster[i]][toCluster[j]]<mini || a==-1) {a=i;b=j;mini=dist[toCluster[i]][toCluster[j]];};
          //Merge two clusters
          nodes[k]=SlotAllocator<Node>::GetSlotAllocator().Allocate();
          new (nodes[k]) Node();
          nodes[k]->setId(k);
          nodes[k]->setAge(mini);
          nodes[k]->setLeft(nodes[toCluster[a]]);
          nodes[k]->setRight(nodes[toCluster[b]]);
          nodes[toCluster[a]]->setFather(nodes[k]);
          nodes[toCluster[b]]->setFather(nodes[k]);
          for (unsigned int i=0;i<toCluster.size();i++) {
            dist[toCluster[i]][k]=(dist[toCluster[i]][toCluster[a]]*toClusterSize[a]+dist[toCluster[i]][toCluster[b]]*toClusterSize[b])/(toClusterSize[a]+toClusterSize[b]);
            dist[k][toCluster[i]]=dist[toCluster[i]][k];
          }
          toClusterSize[a]=toClusterSize[a]+toClusterSize[b];
          toClusterSize[b]=toClusterSize.back();
          toClusterSize.pop_back();
          toCluster[a]=k;
          toCluster[b]=toCluster.back();
          toCluster.pop_back();
        k++;
        }

      //Rescale tree to fit coalescent expectation
      root=nodes[n+n-2];
      double curage=root->getAge();
      double expectedage=2.0*(1.0-1.0/n);
      for (int i=0;i<n+n-1;i++) nodes[i]->setAge(nodes[i]->getAge()*expectedage/curage);
      computeTTotal();
    }

  } // end namespace weakarg
