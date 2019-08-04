#include<iostream>


class TreeNode {
    private:
        TreeNode* parent;
        char value;

    public: 
        TreeNode(char value) {
            this->value = value;
        }

        TreeNode(TreeNode* parent, char value) {
            this->parent = parent;
            this->value = value;   
        }

        char getValue() {
            return this->value;
        }

        TreeNode* getParent() {
            return this->parent;
        }

};

int main(int argc, char** argv) {
    TreeNode t = TreeNode('c');
    std::cout << t.getValue();
}
