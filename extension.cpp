/**
  * File:            extension.cpp
  *
  * Section Leader:  Jonathan Kula
  * Author:          Edouard Des Parois Perrault
  * Date:            Summer 2020
  * Course:          CS106B
  * Quarter:         Summer
  *
  * Summary of File:
  *
  *     This file is an implementation of Shannon entropy to prove that Huffman
  *     encoding produces the shortest possible encoding for a given string. It also
  *     includes some explanations of Shannon entropy. Shannon Entropy is the result of
  *     a formula created by Claude Shannon at Bell Labs inc. in the US. Note that this
  *     file makes use of some functions found in the huffman.cpp file, hence the function
  *     prototypes in huffman.h
  *
  */

#include "extension.h"
#include "huffman.h"
#include "map.h"
#include "vector.h"
#include "testing/SimpleTest.h"
#include <math.h>
using namespace std;
typedef Map<char, double> ProbTable;

ProbTable buildProbTable(FreqTable ft, int messageLength);
double shannonForumula(double pi);
/**
 * This function takes in a string and returns the Shannon entropy for
 * that string. The Shannon entropy is the minimum amount of bits needed
 * per symbol within that string. It can be used to prove that a given Huffman
 * tree makes use of the smallest possible encoding per symbol.
 *
 * This function assumes ft and pt are empty.
 *
 * This function returns a decimal, though keep in mind that, in practice,
 * it is impossible for a bit to be divided into two.
 */
double entropy(string message, FreqTable& ft, ProbTable& pt) {
    // Build a frequency table
    ft = buildFrequencyTable(message);
    // Map out the probabilities
     pt = buildProbTable(ft, message.length());
    double totalEntropy = 0.0;
    for (char ch : pt) {
        totalEntropy += shannonForumula(pt[ch]);
    }
    return totalEntropy;
}

/**
 * This function takes in a string as well as a Huffman tree, and
 * determines whether or not this tree creates codes that use the smallest possible
 * amount of bits per symbol. This function makes use of the formula behind
 * Shannon entropy. In theory, all messages should return true (so long as the corresponding tree is provided),
 * as this is more of an informal proof than a test. That said, because of the fact that Shannon entropy
 * is more precise than the actual length of the Huffman string, it may result in
 * false negatives. This function is not destructive and will not destroy the tree (the tree must be
 * explicitly deallocated afterwards).
 */

ProbTable buildProbTable(FreqTable ft, int messageLength);
double shannonForumula(double pi);
typedef Map<char, int> CodeTable;
CodeTable createCodeTable(EncodingTable et);
double calculateHuffLength(CodeTable& ct, FreqTable& ft);
bool check1(double messageEntropy, CodeTable ct, ProbTable pt);
bool check2(string message, double messageEntropy, CodeTable ct, FreqTable ft);

bool isSmallestPossibleEncoding(string message, EncodingTreeNode* tree) {
    // Declate the necessary variables
    EncodingTable et = createEncodingTable(tree);
    CodeTable ct = createCodeTable(et);
    FreqTable ft;
    ProbTable pt;
    double messageEntropy = entropy(message, ft, pt);
    // Complete both checks
    return check1(messageEntropy, ct, pt) && check2(message, messageEntropy, ct, ft);
}
ProbTable buildProbTable(FreqTable ft, int messageLength) {
    ProbTable pt;
    for (char ch : ft) {
        // 0 is the
        double pi = ft[ch];
        pt[ch] = pi/messageLength;
    }
    return pt;
}
double shannonForumula(double pi) {
    return pi * log2(1 / pi);
}
CodeTable createCodeTable(EncodingTable et) {
    CodeTable ct;
    for (char ch : et) {
        ct[ch] = et[ch].size();
    }
    return ct;
}
double calculateHuffLength(CodeTable& ct, FreqTable& ft) {
    double huffmanSize = 0;
    for (char ch : ct) {
        // Multiply the length of the code by the amount of times is appears
        // inside of the string.
        huffmanSize += ct[ch] * ft[ch];
    }
    return huffmanSize;
}
bool check1(double messageEntropy, CodeTable ct, ProbTable pt) {
    // Multiply the probability of each character by the length of its encoding and add up all characters.
    // This should be equal to the shannon entropy of the text.
    double treeEntropy = 0;
    for (char ch : ct) {
        treeEntropy += ct[ch] * pt[ch];
    }
    // The ceilling function is used to adjust for the fact that
    // there is more precision in Shannon entropy than the actual length
    // generated via Huffman encoding, as bits cannot be split into two.
    return ceil(treeEntropy) == ceil(messageEntropy);
}

bool check2(string message, double messageEntropy, CodeTable ct, FreqTable ft) {
    const int ERROR_MARGIN = 5;
    // We can also solve for the minmum possible encoding length. If there is a minimum of x bits per symbol, and
    // there are y symbols, then Shannon Entropy tells us that the minimum possible encoding length is xy. If the
    // encoded messahe is not xy, there is a problem.
    double minimumPossibleLength = messageEntropy * message.length(); // Shannon entropy tells us that this is the minimum possible length
    double huffLength = calculateHuffLength(ct, ft);
    // Note that the two do not line up. This is because of the fact that the Shannon entropy can be a double
    // of significant precision, while the actual lengtgh of the message cannot, as it is not possible for there to be only a portion of
    // a bit. In other words, Shannon entropy is the theoretically smallest encoding size, though this encoding size is not realistically achievable.
    // Consequently, both numbers do not actually line up exatcly. To correct this, I have put a margin of error of 5. That said, I
    // suppose, as the size of the string increases, the error becomes more and more significant. This margin may need to be adjusted depending
    // on the message size.
    int lowerBound = huffLength - ERROR_MARGIN;
    // Shannon entropy will yield a smaller number than the actual encoding. It may be helpful to think of the actual encoding as rounded up.
    return lowerBound <= minimumPossibleLength && minimumPossibleLength <= huffLength;
}

STUDENT_TEST("Get the entropy for a given string") {
    // Define Variables
    FreqTable ft;
    ProbTable pt;
    // A few examples
    EXPECT_EQUAL(entropy("AAAABBCD", ft, pt),1.75);
    // The following are very precies and are thus rounded up for simplicity
    EXPECT_EQUAL(ceil(entropy("Know then thyself, presume not god to scan", ft, pt)), 5); // Alexander Pope
    EXPECT_EQUAL(ceil(entropy("To be or not to be, that is the question", ft, pt)), 4); // Shakespeare
}

STUDENT_TEST("testing the first check") {
    string message = "AAAABBCD";
    EncodingTreeNode* tree = buildHuffmanTree(message);
    EncodingTable et = createEncodingTable(tree);
    CodeTable ct = createCodeTable(et);
    FreqTable ft;
    ProbTable pt;
    double messageEntropy = entropy(message, ft, pt);
    EXPECT(check1(messageEntropy, ct, pt));
    deallocateTree(tree);
}

STUDENT_TEST("testing the second check") {
    string message = "AAAABBCD";
    EncodingTreeNode* tree = buildHuffmanTree(message);
    EncodingTable et = createEncodingTable(tree);
    CodeTable ct = createCodeTable(et);
    FreqTable ft;
    ProbTable pt;
    double messageEntropy = entropy(message, ft, pt);
    EXPECT(check2(message, messageEntropy, ct, ft));
    deallocateTree(tree);
}

STUDENT_TEST("Testing the entropies of a few strings") {
    string message = "AAAABBCD";
    EncodingTreeNode* tree = buildHuffmanTree(message);
    EXPECT(isSmallestPossibleEncoding(message, tree));
    deallocateTree(tree);
    message = "What reason weaves, by passion is undone"; // More Alexander Pope
    tree = buildHuffmanTree(message);
    EXPECT(isSmallestPossibleEncoding(message, tree));
    deallocateTree(tree);
}

STUDENT_TEST("Testing the entropy using an incorrect tree") {
    // Create an invalid tree
    // (the tree that was created has characters that are not in the string)
    // In some cases, a tree can contain characters thar are not in the string and still
    // be optimal. This is not the case here, however.
    EncodingTreeNode* leaf1 = new EncodingTreeNode('O');
    EncodingTreeNode* leaf2 = new EncodingTreeNode('A');
    EncodingTreeNode* leaf3 = new EncodingTreeNode('C');
    EncodingTreeNode* leaf4 = new EncodingTreeNode('D');
    EncodingTreeNode* internalNode1 = new EncodingTreeNode(leaf1, leaf2);
    EncodingTreeNode* internalNode2 = new EncodingTreeNode(leaf3, leaf4);
    EncodingTreeNode* tree = new EncodingTreeNode(internalNode1, internalNode2);
    // Test the message out
    string message = "OAO";
    EXPECT(!isSmallestPossibleEncoding(message, tree));
    deallocateTree(tree);
}
