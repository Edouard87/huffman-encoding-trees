#pragma once

#include "bits.h"
#include "treenode.h"
#include "queue.h"
#include <string>
#include "priorityqueue.h"
#include "map.h"
using namespace std;

// Helper functions, structs, and typedefs

typedef Queue<Bit> BitStream;
string decodeTextHelper(EncodingTreeNode* tree, Queue<Bit>& messageBits);
typedef Map<char, int> FreqTable;
typedef PriorityQueue<EncodingTreeNode*> Forest;
FreqTable buildFrequencyTable(string text);
Forest buildForest(FreqTable& ft);
EncodingTreeNode* processForest(Forest f); // Turn the Forest into a tree
struct SubTree;
struct TreePair;
typedef Map<char, Vector<Bit>> EncodingTable;
void createEncodingTableHelper(EncodingTreeNode* tree, Vector<Bit>& soFar,  EncodingTable& et);
EncodingTable createEncodingTable(EncodingTreeNode* tree);

// Required prototypes
// Your function implementations must match these without changes

void deallocateTree(EncodingTreeNode* t);
EncodingTreeNode* buildHuffmanTree(std::string messageText);

std::string decodeText(EncodingTreeNode* tree, Queue<Bit>& messageBits);
Queue<Bit> encodeText(EncodingTreeNode* tree, std::string messageText);

void flattenTree(EncodingTreeNode* tree, Queue<Bit>& treeBits, Queue<char>& treeLeaves);
EncodingTreeNode* unflattenTree(Queue<Bit>& treeBits, Queue<char>& treeLeaves);

EncodedData compress(std::string messageText);
std::string decompress(EncodedData& data);
