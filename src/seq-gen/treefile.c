/*  
   Sequence Generator - seq-gen, version 1.3.2
   Copyright (c)1996-2004, Andrew Rambaut & Nick Grassly
   Department of Zoology, University of Oxford			
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://evolve.zoo.ox.ac.uk/software/Seq-Gen/
   email: andrew.rambaut@zoo.ox.ac.uk
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "global.h"
#include "treefile.h"

char treeErrorMsg[256];
int treeError;

TNode *avail=NULL;
long usedAvail=0;
long usedMalloc=0;

/* prototypes */

TNode *NewNode(TTree *tree);
void InitTree(TTree *tree);
void CheckCapacity(TTree *tree, int required);
void DisposeNode(TNode *aNode);
void DisposeTreeNodes(TTree *tree);
void FreeNodes(void);

char ReadToNextChar(FILE *fv);
void ReadUntil(FILE *fv, char stopChar, char *what);
TNode *ReadTip(FILE *fv, char ch, TTree *tree, int numNames, char **names);
TNode *ReadNode(FILE *fv, TTree *tree, int numNames, char **names, int detectPolytomies);
TNode *ReadBranch(FILE *fv, TTree *tree, int numNames, char **names);
void WriteNode(FILE *fv, TTree *tree, TNode *node);


/* functions */

/*----------*/
TNode *NewNode(TTree *tree)
{
	TNode *node;

	if ( avail!=NULL ) {
		node=avail;
		avail=avail->next;
		
		usedAvail++;
	} else {
		if ( (node=malloc(sizeof(TNode)))==NULL ) {
			strcpy(treeErrorMsg, "Out of memory");
			return NULL;
		}
		usedMalloc++;
	}
	node->branch0=node->branch1=node->branch2=NULL;
	node->length0=node->length1=node->length2=0.0;
	node->param=0.0;
	node->tipNo=-1;

	node->sequence = NULL;
	
	node->next=tree->nodeList;	
	tree->nodeList=node;
	
	tree->numNodes++;

	return node;
} /* NewNode */


/*----------*/
void InitTree(TTree *tree)
{	
	tree->root=NULL;
	tree->nodeList=NULL;
	tree->numNodes=0;
	tree->numTips=0;
	tree->totalLength=0.0;
	tree->rooted=0;
	tree->lengths=-1;
} /* InitTree */


/*----------*/
TTree *NewTree(int max_tips)
{
	TTree *tree;
	
	if ( (tree=(TTree *)malloc(sizeof(TTree)))==NULL ) {
		strcpy(treeErrorMsg, "Out of memory creating tree.");
		return NULL;
	}
	tree->capacity=0;
	tree->names = NULL;	//WCC
	tree->tips = NULL;	//WCC
//WCC	CheckCapacity(tree, 1000);
	CheckCapacity(tree, max_tips);
	
	InitTree(tree);
	
	return tree;
} /* NewTree */


/*----------*/
void CheckCapacity(TTree *tree, int required)
{
	int i;
	int newCapacity = tree->capacity;
	char **newNames;
	TNode **newTips;
	
	while (newCapacity < required) {
		newCapacity += 1000;
	}
	
	newNames = (char**)CAllocMem(sizeof(char*)*newCapacity, "newNames", "CheckCapacity", 0);
	newTips = (TNode**)CAllocMem(sizeof(TNode*)*newCapacity, "newTips", "CheckCapacity", 0);
	
	for (i = 0; i < tree->capacity; i++) {
		newNames[i] = tree->names[i];
		newTips[i] = tree->tips[i];
	}
	
	for (i = tree->capacity; i < newCapacity; i++) {
		newNames[i] = NULL;
		newTips[i] = NULL;
	}
	
//WCC	free(tree->names);
//WCC	free(tree->tips);
	if(tree->names != NULL && tree->tips != NULL){
		free(tree->names);
		free(tree->tips);
	}
	
	tree->names = newNames;
	tree->tips = newTips;
	tree->capacity = newCapacity;
} /* CheckCapacity */


/*----------*/
void DisposeNode(TNode *aNode)
{
	aNode->next=avail;
	avail=aNode;
} /* DisposeNode */
			

/*----------*/
void DisposeTreeNodes(TTree *tree)
{
	TNode *P, *O;
	
	if ( tree ) {
		P=tree->nodeList;
		while (P!=NULL) {
			O=P;
			P=P->next;
			DisposeNode(O);
		}
		tree->nodeList=NULL;
	}
} /* DisposeTreeNodes */


/*----------*/
void DisposeTree(TTree *tree)
{
	if (tree) {
		DisposeTreeNodes(tree);
		InitTree(tree);
	}
} /* DisposeTree */


/*----------*/
void FreeNodes(void)
{
	TNode *P, *O;
	
	P=avail;
	while (P!=NULL) {
		O=P;
		P=P->next;
		free(O->sequence);	//WCC
		free(O);
	}
} /* FreeNodes */


/*----------*/
void FreeTree(TTree *tree)
{
	int i;	//WCC

	if (tree) {
		DisposeTreeNodes(tree);
		
		/* Add by Wei-Chen Chen. */
		for(i = 0; i < tree->capacity; i++){
			if(tree->names[i] != NULL){
				free(tree->names[i]);
			}
		}
		free(tree->names);
		free(tree->tips);

		free(tree);
	}
	FreeNodes();
} /* FreeTree */


/*----------*/
//WCC void WriteAvailInfo()
void WriteAvailInfo(void)
{
	TNode *P;
	int count;
	
	count=0;
	P=avail;
	while (P!=NULL) {
		P=P->next;
		count++;	
	}

//WCC	fprintf(stderr, "Avail: %d nodes - availed: %ld, malloced: %ld\n", count, usedAvail, usedMalloc);
	REprintf("Avail: %d nodes - availed: %ld, malloced: %ld\n", count, usedAvail, usedMalloc);
} /* WriteAvailInfo */


int CountTrees(FILE *fv)
{
	int n;
		
	if (fv==NULL)
		return 0;
		
	n=0;
	while (!feof(fv)) {
		if (fgetc(fv)==';')
			n++;
	}
	rewind(fv);
	
	return n;
}

char ReadToNextChar(FILE *fv)
{
	char ch;
	
	ch=fgetc(fv);
	while (!feof(fv) && isspace(ch)) 
		ch=fgetc(fv);

	return ch;
}


void ReadUntil(FILE *fv, char stopChar, char *what)
{
	char ch;
	
	ch=fgetc(fv);
	while (!feof(fv) && ch!=stopChar && 
			ch!='(' && ch!=',' && ch!=':' && ch!=')' && ch!=';') 
		ch=fgetc(fv);

	if (feof(fv) || ch!=stopChar) {
		/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
		// sprintf(treeErrorMsg, "%s missing", what);
		snprintf(treeErrorMsg, 256, "%s missing", what);
		treeError=1;
	}
}


TNode *ReadTip(FILE *fv, char ch, TTree *tree, int numNames, char **names)
{
	int i;
	char *P;
	char name[256];
	
	TNode *node;
	
	node=NewNode(tree);
	i=0;
	
	P=name;
	while (!feof(fv) && ch!=':' && ch!=',' && ch!=')' && i<MAX_NAME_LEN) {
		if (!isspace(ch)) {
			*P=ch;
			i++;P++;
		}
		ch=fgetc(fv);
	}
	*P='\0';

	CheckCapacity(tree, tree->numTips);
	
	if (numNames == 0) {
		node->tipNo=tree->numTips;
		if (tree->names[node->tipNo]==NULL) {
			if ( (tree->names[node->tipNo]=(char *)malloc(MAX_NAME_LEN+1))==NULL ) {
				strcpy(treeErrorMsg, "Out of memory creating name.");
				return NULL;
			}
		}
		strcpy(tree->names[node->tipNo], name);
	} else {
	/* we already have some names so just look it up...*/
		i = 0;
		while (i < numNames && strcmp(name, names[i]) != 0)
			i++;

		if (i == numNames) {
			/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
			// sprintf(treeErrorMsg, "Taxon names in trees for different partitions do not match.");
			snprintf(treeErrorMsg, 256, "Taxon names in trees for different partitions do not match.");
			return NULL;
		}
		
		node->tipNo=i;
	}
	
	tree->tips[node->tipNo]=node;
	tree->numTips++;

	while (!feof(fv) && ch!=':' && ch!=',' && ch!=')') 
		ch=fgetc(fv);

	if (feof(fv)) {
		/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
		// sprintf(treeErrorMsg, "Unexpected end of file");
		snprintf(treeErrorMsg, 256, "Unexpected end of file");
		return NULL;
	}
	ungetc(ch, fv);
	
	return node;
}


TNode *ReadNode(FILE *fv, TTree *tree, int numNames, char **names, int detectPolytomies)
{
	TNode *node, *node2;
	char ch;
	
	if ((node=NewNode(tree))==NULL)
		return NULL;

	if ((node2=ReadBranch(fv, tree, numNames, names))==NULL)
		return NULL;
	node->branch1=node2;
	node2->branch0=node;
	node->length1=node2->length0;
	ReadUntil(fv, ',', "Comma");
	if (treeError)
		return NULL;
	
	if ((node2=ReadBranch(fv, tree, numNames, names))==NULL)
		return NULL;
	node->branch2=node2;
	node2->branch0=node;
	node->length2=node2->length0;
	
	ch=fgetc(fv);
	while (!feof(fv) && ch!=':' && ch!=',' && ch!=')' && ch!=';') 
		ch=fgetc(fv);
		
	if (detectPolytomies && ch==',') {
/*WCC
		fprintf(stderr, "This tree contains nodes which aren't bifurcations. Resolve the node\n");
		fprintf(stderr, "with zero branch lengths to obtain correct results. This can be done\n");
		fprintf(stderr, "with a program called TreeEdit: http://evolve.zoo.ox.ac.uk/software/TreeEdit\n");
*/
		REprintf("This tree contains nodes which aren't bifurcations. Resolve the node\n");
		REprintf("with zero branch lengths to obtain correct results. This can be done\n");
		REprintf("with a program called TreeEdit: http://evolve.zoo.ox.ac.uk/software/TreeEdit\n");
		exit(0);
	}

	if (feof(fv)) {
		/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
		// sprintf(treeErrorMsg, "Unexpected end of file");
		snprintf(treeErrorMsg, 256, "Unexpected end of file");
		return NULL;
	}
	ungetc(ch, fv);
	
	return node;
}


TNode *ReadBranch(FILE *fv, TTree *tree, int numNames, char **names)
{
	char ch;
	double len, param=0.0;
	TNode *node;
	
	ch=ReadToNextChar(fv);
	if (ch=='(') {	/* is a node */
		node=ReadNode(fv, tree, numNames, names, 1);
		ReadUntil(fv, ')', "Closing bracket");
		if (treeError)
			return NULL;
	} else {		/* is a tip */
		node=ReadTip(fv, ch, tree, numNames, names);
	}
	
	ch=ReadToNextChar(fv);
	if (ch==':') {
		if (tree->lengths==0) {
			/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
			// sprintf(treeErrorMsg, "Some branches don't have branch lengths");
			snprintf(treeErrorMsg, 256, "Some branches don't have branch lengths");
			return NULL;
		} else 
			tree->lengths=1;
			
		if (fscanf(fv, "%lf", &len)!=1) {
			/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
			// sprintf(treeErrorMsg, "Unable to read branch length");
			snprintf(treeErrorMsg, 256, "Unable to read branch length");
			return NULL;
		}

		ch=ReadToNextChar(fv);
		if (ch=='[') {
			if (fscanf(fv, "%lf", &param)!=1) {
				/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
				// sprintf(treeErrorMsg, "Unable to read branch parameter");
				snprintf(treeErrorMsg, 256, "Unable to read branch parameter");
				return NULL;
			}
			ReadUntil(fv, ']', "Close square bracket");
		} else
			ungetc(ch, fv);
	} else {
		if (tree->lengths==1) {
			/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
			// sprintf(treeErrorMsg, "Some branches don't have branch lengths");
			snprintf(treeErrorMsg, 256, "Some branches don't have branch lengths");
			return NULL;
		} else 
			tree->lengths=0;
	
		len=0.0;
		ungetc(ch, fv);
	}
	node->length0=len;
	node->param=param;
	
	return node;
}	


void ReadTree(FILE *fv, TTree *tree, int treeNum, int numNames, char **names, 
				int *outNumSites, double *outRelRate)
{
	char ch;
	TNode *P;
		
	treeError=0;
	tree->numNodes=0;
	tree->numTips=0;
	tree->rooted=1;
	tree->lengths=-1;
	
	(*outRelRate) = 1.0;
	
	ch=fgetc(fv);
	while (!feof(fv) && ch!='(' && ch!='[') 
		ch=fgetc(fv);

	if (ch == '[') {
		if (fscanf(fv, "%d", outNumSites)!=1) {
			/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
			// sprintf(treeErrorMsg, "Unable to read partition length");
			snprintf(treeErrorMsg, 256, "Unable to read partition length");
			exit(0);
		}

		ch=fgetc(fv);
		while (!feof(fv) && ch!=',' && ch!='(') 
			ch=fgetc(fv);
			
		if (ch == ',') {
			if (fscanf(fv, "%lf", outRelRate)!=1) {
				/* R-devel on around Dec. 24, 2022 starting to warn the line below. */
				// sprintf(treeErrorMsg, "Unable to read partition relative rate");
				snprintf(treeErrorMsg, 256, "Unable to read partition relative rate");
				exit(0);
			}

			ch=fgetc(fv);
			while (!feof(fv) && ch!='(') 
				ch=fgetc(fv);
		}
	}
	
	if (ch!='(' || (tree->root=ReadNode(fv, tree, numNames, names, 0))==NULL) {
//WCC		fprintf(stderr, "Error reading tree number %d: %s.\n", treeNum, treeErrorMsg);
		REprintf("Error reading tree number %d: %s.\n", treeNum, treeErrorMsg);
		exit(0);
	}
	
	ch=fgetc(fv);
	while (!feof(fv) && ch!=',' && ch!=')' && ch!=';') 
		ch=fgetc(fv);

	if (ch==',') {		/* the tree is unrooted */
		tree->rooted=0;
		if ((tree->root->branch0=ReadBranch(fv, tree, numNames, names))==NULL) {
//WCC			fprintf(stderr, "Error reading tree number %d: %s.\n", treeNum, treeErrorMsg);
			REprintf("Error reading tree number %d: %s.\n", treeNum, treeErrorMsg);
			exit(0);
		}
		tree->root->branch0->branch0=tree->root;
		tree->root->length0=tree->root->branch0->length0;
	}
	
	tree->totalLength=0.0;
	
	if (tree->rooted) {
		P=tree->root;
		while (P!=NULL) {
			tree->totalLength+=P->length0;
			P=P->branch1;
		}
	}
}


int IsTreeAvail(FILE *fv)
{
	char ch;

	ch=fgetc(fv);
	while (!feof(fv) && ch!='(' && ch!='[') 
		ch=fgetc(fv);
		
	if (ch=='(' || ch=='[')
		ungetc(ch, fv);
		
	return (!feof(fv));
}


void WriteNode(FILE *fv, TTree *tree, TNode *node)
{
	if (node->tipNo==-1) {
		fputc('(', fv);				
		WriteNode(fv, tree, node->branch1);
		fputc(',', fv);			
		WriteNode(fv, tree, node->branch2);
		fputc(')', fv);				
	} else {
		fprintf(fv, "%s", tree->names[node->tipNo]);
	}
	if (tree->lengths)
		fprintf(fv, ":%.6f", node->length0);
}


void WriteTree(FILE *fv, TTree *tree)
{
	fputc('(', fv);	
	if (tree->rooted) {
		WriteNode(fv, tree, tree->root->branch1);
		fputc(',', fv);				
		WriteNode(fv, tree, tree->root->branch2);
	} else {
		WriteNode(fv, tree, tree->root->branch1);
		fputc(',', fv);				
		WriteNode(fv, tree, tree->root->branch2);
		fputc(',', fv);	
		WriteNode(fv, tree, tree->root->branch0);
	}
	fprintf(fv, ");\n");
}


void UnrootRTree(TTree *tree) 
/* Used to unroot a rooted tree */
{
//WCC	TNode *P, *Q, *R, *newRoot;
	TNode *P, *Q, *R = NULL, *newRoot;
	double len, len2;
		
	if (!tree->rooted || tree->numTips<3)
		return;
	
	P=tree->tips[0];
	Q=P->branch0; newRoot=Q;
	while (Q!=tree->root) {
		R=Q->branch0;
		len=Q->length0;
		if (Q->branch1==P) {
			len2=Q->length1;
			Q->branch1=R;
			Q->length1=len;
		} else {
			len2=Q->length2;
			Q->branch2=R;
			Q->length2=len;
		}
		Q->branch0=P;
		Q->length0=len2;

		P=Q;
		Q=R;
	}

	len=R->length1+R->length2;
	if (R->branch1==P)
		Q=R->branch2;
	else
		Q=R->branch1;

	Q->branch0=P;
	Q->length0=len;
	
	if (P->branch1==R) {
		P->branch1=Q;
		P->length1=len;
	} else {
		P->branch2=Q;
		P->length2=len;
	}
	
	tree->root=newRoot;
	
	DisposeNode(R);
	
	tree->rooted=0;
}


void RerootUTree(TTree *tree, int tip) 
/* Used to reroot an unrooted tree */
/* This may sound strange but all it does is move
   the root trichotomy  */
{
	TNode *P, *Q, *R, *newRoot;
	double len, len2;
		
	if (tree->rooted)
		return;
	
	P=tree->tips[tip];
	Q=P->branch0; newRoot=Q;
	while (P!=tree->root) {
		R=Q->branch0;
		len=Q->length0;
		if (Q->branch1==P) {
			len2=Q->length1;
			Q->branch1=R;
			Q->length1=len;
		} else {
			len2=Q->length2;
			Q->branch2=R;
			Q->length2=len;
		}
		Q->branch0=P;
		Q->length0=len2;

		P=Q;
		Q=R;
	}
	 	
	tree->root=newRoot;
}


