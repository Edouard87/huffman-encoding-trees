/**
  * File:            huffman.cpp
  *
  * Section Leader:  Jonathan Kula
  * Author:          Edouard Des Parois Perrault
  * Date:            Summer 2020
  * Course:          CS106B
  * Quarter:         Summer
  *
  * Summary of File:
  *
  *     This file contains my implementation of the huffman coding algorithm. It makes use
  *     of various provided functions as well as additional helper functions and is based on
  *     what was taught in lecture.
  *
  *     So that some of the functions in this file could be implemented in my extension, I have
  *     added my helper function prototypes and my structs to huffman.h. Please ensure that the huffman.h
  *     header file I have submitted is being used.
  *
  */

#include "bits.h"
#include "treenode.h"
#include "huffman.h"
#include "map.h"
#include "vector.h"
#include "priorityqueue.h"
#include "strlib.h"
#include "testing/SimpleTest.h"
using namespace std;

/**
 * Given a Queue<Bit> containing the compressed message bits and the encoding tree
 * used to encode those bits, decode the bits back to the original message text.
 *
 * You can assume that tree is a well-formed non-empty encoding tree and
 * bits queue contains a valid sequence of encoded bits.
 *
 * The first parameter is the root node of the encoding tree and the second parameter is
 * a queue of the message to be encoded using that tree. The tree should be a Huffman tree.
 *
 * In addition, note that this function is capable of recognizing certain incorrect
 * inputs. If it attempts to decode and the supplied encoding is too long, it will throw
 * an error. Similarily, if it expects a longer tree, it will throw an error.
 */

string decodeText(EncodingTreeNode* tree, Queue<Bit>& messageBits) {
    string decodedText;
    while (!messageBits.isEmpty()) {
        decodedText += decodeTextHelper(tree, messageBits);
    }
    return decodedText;
}

string decodeTextHelper(EncodingTreeNode* tree, Queue<Bit>& messageBits) {
    if (tree->zero == nullptr && tree->one == nullptr) {
        return to_string(tree->ch);
    }
    /*
     * The position of the below line is very important.
     * The leaf node will be returned before it reaches these
     * lines if the input is correct.
     */
    if (messageBits.isEmpty()) {
        error("Invalid Input");
    }
    Bit bit = messageBits.dequeue();
    if (bit == 0) {
        if (tree->zero == nullptr) {
            error("Invalid input");
        }
        tree = tree->zero;
    } else if (bit == 1) {
        if (tree->one == nullptr) {
            error("Invalid input");
        }
        tree = tree->one;
    }
    return decodeTextHelper(tree, messageBits);
}

/**
 * Reconstruct encoding tree from flattened form Queue<Bit> and Queue<char>.
 *
 * You can assume that the input Queues are well-formed and represent
 * a valid encoding tree.
 *
 * The first parameter is a sequence of bits that represents the compressed (flattened) tree. The second
 * parameter is a queue of the characters that are to be placed at the leaf nodes. The result is the root node
 * of the encoding tree that conforms to these two parameters. This function cannot handle malformed inputs.
 */

EncodingTreeNode* unflattenTree(Queue<Bit>& treeBits, Queue<char>& treeLeaves) { // 'pre-order' creation of the tree
    Bit bit = treeBits.dequeue(); // Get the next bit
    if (bit == 0) {
        return new EncodingTreeNode(treeLeaves.dequeue()); // 0 is a leaf node
    } else {
       return new EncodingTreeNode(unflattenTree(treeBits, treeLeaves), unflattenTree(treeBits, treeLeaves)); // 1 is an internal node + two other mystery nodes
    }
}

/**
 * Decompress the given EncodedData and return the original text.
 *
 * You can assume the input data is well-formed and was created by a correct
 * implementation of compress.
 *
 * Your implementation may change the data parameter however you like. There
 * are no requirements about what it should look like after this function
 * returns.
 *
 * Note that this function is destructive, and will empty the queues in the data parameter. It will also
 * automatically deallocate the encoding tree. No additional cleaning must be performed after
 * this function is called.
 */
string decompress(EncodedData& data) {
    EncodingTreeNode* tree = unflattenTree(data.treeBits, data.treeLeaves); // Get the encoding tree
    string decodedtext = decodeText(tree, data.messageBits); // Decode the message using the tree
    deallocateTree(tree); // Deallocate the tree after the message was decoded to avoid memory leaks
    return decodedtext;
}

/**
 * Constructs an optimal Huffman coding tree for the given text, using
 * the algorithm described in lecture.
 *
 * Reports an error if the input text does not contain at least
 * two distinct characters.
 *
 * When assembling larger trees out of smaller ones, make sure to set the first
 * tree dequeued from the queue to be the zero subtree of the new tree and the
 * second tree as the one subtree.
 *
 * This function makes use of several helper functions and structs, defined below. It takes in a
 * single parameter that is the text it is to process. It then creates the optimal Huffman encoding
 * tree for that string and returns its root node. It is not equipped to handle malformed inputs.
 */

EncodingTreeNode* buildHuffmanTree(string text) {
    // Build a frequency table
    FreqTable ft = buildFrequencyTable(text);
    // Make sure there are at least 2 distinct entries
    if (ft.size() <= 1) {
        error("Your string only contains one character. Strings to be encoded must contain at least two distinct characters.");
    }
    // Create a node for each character and enqueue it
    Forest f = buildForest(ft);
    // Process the resulting Forest
    EncodingTreeNode* huffmanTree = processForest(f);
    // Return the completed huffman tree
    return huffmanTree;
}

FreqTable buildFrequencyTable(string text) {
    FreqTable ft;
    for (char ch : text) {
        ft[ch]++; // Incrememt the occurence of that char in the frequency table
    }
    return ft;
}
Forest buildForest(FreqTable& ft) {
    Forest q;
    for (char ch : ft) {
        EncodingTreeNode* leaf = new EncodingTreeNode(ch);
        q.enqueue(leaf, ft[ch]);
    }
    return q;
}
struct SubTree {
    double priority;
    EncodingTreeNode* ptr;
    SubTree() {
        this->priority = 0;
        this->ptr = nullptr;
    }
    SubTree(double priority, EncodingTreeNode* ptr) {
        this->priority = priority;
        this->ptr = ptr;
    }
};
struct TreePair{
    SubTree first;
    SubTree second;
    TreePair(SubTree f, SubTree s) {
        this->first = f;
        this->second = s;
    }
    double sumPriority() {
        return this->first.priority + this->second.priority;
    }
};
EncodingTreeNode* processForest(Forest f) {
    // Not done recursively as I fear this may
    // result in Stack Overflow.
    while (f.size() > 1) { // So long as there is more than one tree in the forest
        TreePair tp = {{f.peekPriority(), f.dequeue()},{f.peekPriority(), f.dequeue()}}; // Bundles both subtrees in one struct
        EncodingTreeNode* internalNode = new EncodingTreeNode(tp.first.ptr, tp.second.ptr); // Creates an internal tree node based on this struct
        f.enqueue(internalNode, tp.sumPriority()); // Puts the internal node back in the queue with a priority that is the sum of both subtrees
    }
    return f.dequeue();
}

/**
 * Given a string and an encoding tree, encode the text using the tree
 * and return a Queue<Bit> of the encoded bit sequence.
 *
 * You can assume tree is a valid non-empty encoding tree and contains an
 * encoding for every character in the text.
 *
 * Uses the tree passed as a first parameter to encode the supplied string of text passed as
 * a second parameter. It returns a queue of bits that represents the encoded text.
 */

Queue<Bit> encodeText(EncodingTreeNode* tree, string text) {
    // A vector here is used because we must combine them
    // Step 1: Create encoding table
    EncodingTable et = createEncodingTable(tree); // note that an encoding table is not the same as a frequency table
    // Step 2: Map each character using the table
    BitStream bits;
    for (int i = 0; i < (int) text.length(); i++) {
        char ch = text.at(i);
        Vector<Bit> code = et[ch]; // Get the code for that charater
        for (Bit bit : code) {
            bits.enqueue(bit);
        }
    }
    return bits;
}

void createEncodingTableHelper(EncodingTreeNode* tree, Vector<Bit>& soFar,  EncodingTable& et) { // pre-order traversal
    // do something
    if (tree->one == nullptr && tree->zero == nullptr) {
        // This is a character
        et[tree->ch] = soFar;
        return;
    }
    // traverse left
    soFar.add(0);
    createEncodingTableHelper(tree->zero, soFar, et);
    soFar.remove(soFar.size() - 1);
    // traverse right
    soFar.add(1);
    createEncodingTableHelper(tree->one, soFar, et);
    soFar.remove(soFar.size() - 1);
}

EncodingTable createEncodingTable(EncodingTreeNode* tree) {
    Vector<Bit> soFar;
    EncodingTable et;
    createEncodingTableHelper(tree, soFar, et);
    return et;
}

/**
 * Flatten the given tree into a Queue<Bit> and Queue<char> in the manner
 * specified in the assignment writeup.
 *
 * You can assume the input Queues are empty on entry to this function.
 *
 * You can assume tree is a valid well-formed encoding tree.
 *
 * This function takes in the root node of an encoding tree as well as a queue of bits that will be populated
 * with bits representing the flattened version of the tree (using a pre-order traversal). The last pramameter will be
 * populated by the characters to be placed in the leaf nodes of the tree, ordered from left to right.
 */
void flattenTree(EncodingTreeNode* tree, Queue<Bit>& treeBits, Queue<char>& treeLeaves) {
    // Do something
    if (tree->one == nullptr && tree->zero == nullptr) {
        // This is a leaf node
        treeBits.enqueue(0);
        treeLeaves.add(tree->ch);
        return;
    }
    // This is an internal node
    treeBits.enqueue(1);
    // Traverse left
    flattenTree(tree->zero, treeBits, treeLeaves);
    // Traverse right
    flattenTree(tree->one, treeBits, treeLeaves);
}

/**
 * Compress the input text using Huffman coding, producing as output
 * an EncodedData containing the encoded message and encoding tree used.
 *
 * Reports an error if the message text does not contain at least
 * two distinct characters.
 *
 * Compresses the input string into a struct of encoded data. This data contains the
 * compressed tree, the array of characters that belong in the tree, and the encoded
 * message using the aforementioned parameters.
 */
EncodedData compress(string messageText) {
    // Create a Huffman tree with the data
    EncodingTreeNode* tree = buildHuffmanTree(messageText); // Handles the not enough variety edge case.
    // Flatten the tree
    BitStream treeBits;
    Queue<char> leaves;
    flattenTree(tree, treeBits, leaves);
    // Encode the message
    BitStream encodedText = encodeText(tree, messageText);
    // Deallocate the tree
    deallocateTree(tree);
    // Bundle everything up and return
    EncodedData data;
    data.treeBits = treeBits;
    data.treeLeaves = leaves;
    data.messageBits = encodedText;
    return data;
}

/* * * * * * Testing Helper Functions Below This Point * * * * * */

EncodingTreeNode* createExampleTree() {
    /* Example encoding tree used in multiple test cases:
     *                *
     *              /   \
     *             T     *
     *                  / \
     *                 *   E
     *                / \
     *               R   S
     */
    /* Variable Declarations */
    EncodingTreeNode* zero;
    EncodingTreeNode* one;
    EncodingTreeNode* parent;
    /* Creating the Tree */
    zero = new EncodingTreeNode('R');
    one = new EncodingTreeNode('S');
    parent = new EncodingTreeNode(zero, one);
    zero = parent;
    one = new EncodingTreeNode('E');
    parent = new EncodingTreeNode(zero, one);
    zero = new EncodingTreeNode('T');
    one = parent;
    return new EncodingTreeNode(zero, one);
}

void deallocateTree(EncodingTreeNode* root) {
    /* Base Case */
    if (root == nullptr) {
        return;
    }
    /* Recursive Case */
    deallocateTree(root->zero);
    deallocateTree(root->one);
    delete root;
}

bool isGood(EncodingTreeNode* a, EncodingTreeNode* b) {
    if (a == nullptr && b == nullptr) {
        return true;
    }
    if (a != nullptr && b != nullptr) {
        return true;
    } else {
        return false;
    }
}

bool areEqual(EncodingTreeNode* a, EncodingTreeNode* b) {
    /* Base Case */
    /* Compare the two characters */
    if (a == nullptr && b == nullptr) {
        return true;
    } else if (a != nullptr && b != nullptr && a->ch == b->ch) {
        return areEqual(a->zero, b->zero) && areEqual(a->one, b->one);
    } else {
        return false;
    }
    // cannot simply compare a to b, as a and b should be different addresses.
    /* Recursive Case */

}

/* * * * * * Test Cases Below This Point * * * * * */

STUDENT_TEST("Create an example tree cand delete it.") {
    EncodingTreeNode* tree = createExampleTree();
    deallocateTree(tree);
}

STUDENT_TEST("Check for equality between two equal trees") {
    EncodingTreeNode* tree = createExampleTree();
    EncodingTreeNode* tree2 = createExampleTree();
    EXPECT(areEqual(tree, tree2));
    deallocateTree(tree);
    deallocateTree(tree2);
}

STUDENT_TEST("Check for equality between two equal subtrees") {
    EncodingTreeNode* tree = createExampleTree();
    EncodingTreeNode* tree2 = createExampleTree();
    EncodingTreeNode* subtree = tree->one->zero;
    EncodingTreeNode* subtree2 = tree->one->zero;
    EXPECT(areEqual(subtree, subtree2));
    deallocateTree(tree);
    deallocateTree(tree2);
}

BitStream stringToBitStream(string input) {
    Vector<int> output1;
    BitStream output2;
    for (int i = 0; i < (int) input.size(); i++) {
        int currentNum = input.at(i) - 48;
        output1.add(currentNum);
    }
    for (int num : output1) {
        output2.enqueue(num);
    }
    return output2;
}

STUDENT_TEST("Check for equality between two different subtrees") {
    EncodingTreeNode* tree = createExampleTree();
    EncodingTreeNode* tree2 = createExampleTree();
    EncodingTreeNode* subtree = tree->one;
    EncodingTreeNode* subtree2 = tree->one->zero->one;
    EXPECT(!areEqual(subtree, subtree2));
    deallocateTree(tree);
    deallocateTree(tree2);
}

STUDENT_TEST("Test the string to bitstream helper") {
    BitStream bitstream = stringToBitStream("0100101111");
    EXPECT_EQUAL(bitstream.dequeue(), 0);
    EXPECT_EQUAL(bitstream.dequeue(), 1);
}

STUDENT_TEST("Test the decode function") {
    BitStream bitstream = stringToBitStream("010010111");
    EncodingTreeNode* tree = createExampleTree();
    string decodedText = decodeText(tree, bitstream);
    EXPECT_EQUAL(decodedText, "TRSE");
    deallocateTree(tree);
}

STUDENT_TEST("The decode function can handle incorrect inputs") {
    BitStream bitstream = stringToBitStream("0100101111");
    EncodingTreeNode* tree = createExampleTree();
    EXPECT_ERROR(decodeText(tree, bitstream));
    deallocateTree(tree);
}

STUDENT_TEST("Unflatten a falttened tree") {
    BitStream treeBits = stringToBitStream("1101000"); // Example tree
    Queue<char> leaves = {'N','M','S','O'};
    EncodingTreeNode* tree = unflattenTree(treeBits, leaves);
    deallocateTree(tree);
}

STUDENT_TEST("Unflatten another falttened tree") {
    BitStream treeBits = stringToBitStream("1100100"); // My Tree
    Queue<char> leaves = {'A','B','C','D'};
    EncodingTreeNode* tree = unflattenTree(treeBits, leaves);
    deallocateTree(tree);
}

STUDENT_TEST("Decompress") {
    EncodedData data;
    data.treeBits = stringToBitStream("1100100"); // My Tree
    data.treeLeaves = {'A','B','C','D'};
    data.messageBits = stringToBitStream("00011011");
    EXPECT_EQUAL(decompress(data), "ABCD");
}

STUDENT_TEST("Creating an encoding table") {
    BitStream treeBits = stringToBitStream("1100100");
    Queue<char> leaves = {'A','B','C','D'};
    EncodingTreeNode* tree = unflattenTree(treeBits, leaves);
    EncodingTable et = createEncodingTable(tree);
    Vector<Bit> ACode = {0, 0},
            BCode = {0, 1},
            CCode = {1, 0},
            DCode = {1, 1};
    EXPECT_EQUAL(et['A'], ACode);
    EXPECT_EQUAL(et['B'], BCode);
    EXPECT_EQUAL(et['C'], CCode);
    EXPECT_EQUAL(et['D'], DCode);
    deallocateTree(tree);
}

STUDENT_TEST("Encoding a message") {
    string message = "ABCD";
    BitStream treeBits = stringToBitStream("1100100");
    Queue<char> leaves = {'A','B','C','D'};
    EncodingTreeNode* tree = unflattenTree(treeBits, leaves);
    BitStream bits = encodeText(tree, message);
    EXPECT_EQUAL(bits.dequeue(),0);
    EXPECT_EQUAL(bits.dequeue(),0);
    deallocateTree(tree);
}

STUDENT_TEST("Flatten a tree") {
    // Create a tree
    string message = "ABCD";
    BitStream treeBits = stringToBitStream("1100100");
    BitStream oldTreeBits = treeBits; // Make a copy; the flattening is destructive
    Queue<char> leaves = {'A','B','C','D'};
    Queue<char> oldLeaves = leaves; // Make a copy; the flattening is destructive
    EncodingTreeNode* tree = unflattenTree(treeBits, leaves);
    // Flatten the tree
    BitStream newTreeBits;
    Queue<char> newLeaves;
    flattenTree(tree, newTreeBits, newLeaves);
    // Compare the results
    EXPECT_EQUAL(newLeaves, oldLeaves);
    EXPECT_EQUAL(newTreeBits, oldTreeBits);
    // Deallocate everything else
    deallocateTree(tree);
}

STUDENT_TEST("Build frequency table") {
    string text = "weclome to cs106B";
    FreqTable ft = buildFrequencyTable(text);
    EXPECT_EQUAL(ft['w'],1);
    EXPECT_EQUAL(ft['e'],2);
    EXPECT_EQUAL(ft['1'],1);
    text = "know then thyself"; // Quote by Alexander Pope
    ft = {};
    ft = buildFrequencyTable(text);
    EXPECT_EQUAL(ft['n'],2);
    EXPECT_EQUAL(ft['e'],2);
}

STUDENT_TEST("Create a Forest") {
    string text = "know then thyself"; // Quote by Alexander Pope
    FreqTable ft = buildFrequencyTable(text);
    Forest q = buildForest(ft);
    EXPECT_EQUAL(q.peekPriority(), 1);
    while (!q.isEmpty()) { // Deallocate
        delete q.dequeue();
    }
}

STUDENT_TEST("Turn the forest into a Huffman Tree") {
    string text = "heeppp";
    FreqTable ft = buildFrequencyTable(text);
    Forest f = buildForest(ft);
    EncodingTreeNode* tree = processForest(f);
    EXPECT_EQUAL(tree->one->ch, 'p'); // It is easier to inspect the full tree in the debugger
    deallocateTree(tree);
}

STUDENT_TEST("Compress a message") {
    string text = "hello";
    EncodedData compressed = compress(text);
    // The following were calculated by hand
    BitStream encodedMessage = stringToBitStream("1111000110");
    BitStream FlattenedTree = stringToBitStream("1010100");
    Queue<char> leaves = {'l','e','o','h'};
    // Compare to those calculated by the computer
    EXPECT_EQUAL(compressed.messageBits,encodedMessage);
    EXPECT_EQUAL(compressed.treeBits,FlattenedTree);
    EXPECT_EQUAL(compressed.treeLeaves,leaves);
}

STUDENT_TEST("Test compress -> decompress on a very large text") {
    /* based on the PROVIDED TEST pertaining on this same subject */
    Vector<string> inputs = {
        "lotsT of$ special^&# characters&",
        "(*) random message."
    };

    for (string input: inputs) {
        EncodedData data = compress(input);
        string output = decompress(data);

        EXPECT_EQUAL(output.size(), input.size());

        /* Don't clobber the output with a huge string if there's a mismatch. */
        EXPECT(input == output);
    }
}

STUDENT_TEST("Throw an error if there is only one character") {
    EXPECT_ERROR(compress("e"));
    EXPECT_ERROR(compress("eeeeeeeee"));
}

/* * * * * Provided Tests Below This Point * * * * */

PROVIDED_TEST("decodeText, small example encoding tree") {
    EncodingTreeNode* tree = createExampleTree(); // see diagram above

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(decodeText(tree, messageBits), "E");

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(decodeText(tree, messageBits), "SET");

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1}; // STREETS
    EXPECT_EQUAL(decodeText(tree, messageBits), "STREETS");

    deallocateTree(tree);
}

PROVIDED_TEST("unflattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  treeBits   = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'T', 'R', 'S', 'E' };
    EncodingTreeNode* tree = unflattenTree(treeBits, treeLeaves);

    EXPECT(areEqual(tree, reference));

    deallocateTree(tree);
    deallocateTree(reference);
}

PROVIDED_TEST("decompress, small example input") {
    EncodedData data = {
        { 1, 0, 1, 1, 0, 0, 0 }, // treeBits
        { 'T', 'R', 'S', 'E' },  // treeLeaves
        { 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1 } // messageBits
    };
    EXPECT_EQUAL(decompress(data), "TRESS");
}

PROVIDED_TEST("encodeText, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    EncodingTreeNode* tree = buildHuffmanTree("STREETTEST");
    EXPECT(areEqual(tree, reference));
    deallocateTree(reference);
    deallocateTree(tree);
}

PROVIDED_TEST("encodeText, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(encodeText(reference, "E"), messageBits);

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(encodeText(reference, "SET"), messageBits);

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1 }; // STREETS
    EXPECT_EQUAL(encodeText(reference, "STREETS"), messageBits);

    deallocateTree(reference);
}

PROVIDED_TEST("flattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  expectedBits   = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> expectedLeaves = { 'T', 'R', 'S', 'E' };

    Queue<Bit>  treeBits;
    Queue<char> treeLeaves;
    flattenTree(reference, treeBits, treeLeaves);

    EXPECT_EQUAL(treeBits,   expectedBits);
    EXPECT_EQUAL(treeLeaves, expectedLeaves);

    deallocateTree(reference);
}

PROVIDED_TEST("compress, small example input") {
    EncodedData data = compress("STREETTEST");
    Queue<Bit>  treeBits    = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeChars   = { 'T', 'R', 'S', 'E' };
    Queue<Bit>  messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0 };

    EXPECT_EQUAL(data.treeBits, treeBits);
    EXPECT_EQUAL(data.treeLeaves, treeChars);
    EXPECT_EQUAL(data.messageBits, messageBits);
}

PROVIDED_TEST("Test end-to-end compress -> decompress") {
    Vector<string> inputs = {
        "HAPPY HIP HOP",
        "The job requires extra pluck and zeal from every young wage earner.",
        ":-) :-D XD <(^_^)>",
    };

    for (string input: inputs) {
        EncodedData data = compress(input);
        string output = decompress(data);

        EXPECT_EQUAL(output.size(), input.size());

        /* Don't clobber the output with a huge string if there's a mismatch. */
        EXPECT(input == output);
    }
}
