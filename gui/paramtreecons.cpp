#include "paramtreecons.h"

ParamTreeCons::ParamTreeCons() {its=0;}

ParamTreeCons::~ParamTreeCons() {}

void ParamTreeCons::account()
{
	its++;
	int n=getTree()->getN();
	for (int i=0;i<n*2-1;i++) {
		Node * node=getTree()->getNode(i);
		QString key(n,'0');
		makeKey(node,&key);
		if (hash.contains(key)) {
			hash[key].d=(hash[key].d*hash[key].i+node->getAge())/(1.0+hash[key].i);
			hash[key].i++;
			}
		else hash[key]=Cluster(1,node->getAge());
	}
}

void ParamTreeCons::makeKey(Node * node,QString * key)
{
	if (node->getLeft()==NULL) key->replace(node->getId(),1,'1');else{
	makeKey(node->getLeft (),key);
	makeKey(node->getRight(),key);
	}
}

void ParamTreeCons::consensus(int cutoff)
{
 vector<QString> toKeep;
 QHashIterator<QString, Cluster> i(hash);
 while (i.hasNext()) {
     i.next();
     if (i.value().i>=1.0*cutoff*its/100) toKeep.push_back(i.key());
 }
 QString key(getTree()->getN(),'1');
 QString res=buildsubtree(key,&toKeep);
 res.append(":0.000000");
 delete(rectree);
 rectree=new RecTree(data->getL(),res.toStdString(),false,false);
}

QString ParamTreeCons::buildsubtree(QString key,vector<QString>* tokeep)
{
int card=key.count('1');
double age=hash[key].d;
vector<QString> subs;
for (int i=card-1;i>0;i--) {
	for (unsigned int j=0;j<tokeep->size();j++) {
		if (tokeep->at(j).count('1')==i) {
			bool ok=true;
			for (int k=0;k<key.length();k++)  if (tokeep->at(j)[k]=='1' && key[k]=='0') {ok=false;break;}
			if (!ok) continue;
			if (i==1) {
				subs.push_back(QString::number(tokeep->at(j).indexOf('1')).append(":").append(QString::number(age)));
			} else {
				subs.push_back(buildsubtree(tokeep->at(j),tokeep).append(":").append(QString::number(max(0.0,age-hash[tokeep->at(j)].d))));
			}
		for (int k=0;k<key.length();k++) if (tokeep->at(j)[k]=='1') key[k]='0';
		}
	}
}
QString res(subs.size()-1,'(');
res.append(subs[0]).append(",");
for (unsigned int i=1;i<subs.size();i++) res.append(subs[i]).append("):0.0,");
res.chop(5);
return res;
}

void ParamTreeCons::consensusExt()
{
 QHash<QString, Cluster> h=hash;
 vector<QString> toKeep;
 while (h.size()>0) {
 QHashIterator<QString, Cluster> i(h);
 i.next();
 QString cur=i.key();
 while (i.hasNext()) {
     i.next();
     if (i.value().i>h[cur].i) cur=i.key();
 }
 bool ok=true;
 for (unsigned int j=0;j<toKeep.size();j++) if (!compat(cur,toKeep[j])) {ok=false;break;};
 if (ok) toKeep.push_back(cur);
 h.remove(cur);
 }
 QString key(getTree()->getN(),'1');
 QString res=buildsubtree(key,&toKeep);
 res.append(":0.000000");
 delete(rectree);
 rectree=new RecTree(data->getL(),res.toStdString(),false,true);
}

bool ParamTreeCons::compat(QString key1,QString key2)
{
	bool ok=true;
	for (int i=0;i<key1.length();i++) if (key1[i]=='1'&&key2[i]=='1') {ok=false;break;}
	if (ok) return true;
	ok=true;
	for (int i=0;i<key1.length();i++) if (key1[i]=='1'&&key2[i]=='0') {ok=false;break;}
	if (ok) return true;
	ok=true;
	for (int i=0;i<key1.length();i++) if (key1[i]=='0'&&key2[i]=='1') {ok=false;break;}
	return ok;
}

