#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <memory>
#include <cassert>

#include "wrap_xml.h"

#include "tinyxml2.h"

namespace mp {

using namespace tinyxml2;

/// This stores either the whole tree,
/// or a subnode.
class TinyXML2Walker
    : public BasicTreeWalker {
  std::unique_ptr<XMLDocument> pdoc_;
  XMLElement* pelem_ {};
public:
  /// Default-construct
  TinyXML2Walker() : pdoc_(std::make_unique<XMLDocument>()) { }

  /// Construct for a node
  TinyXML2Walker(XMLElement* pe) : pelem_(pe) { assert(pe); }

  /// Read file (for the root only)
  bool ReadFile(const char* name) override {
    assert(pdoc_);
    XMLError result = pdoc_->LoadFile(name);
    if (result != XML_SUCCESS) {
      std::cerr << "Error loading XML file: " << name << std::endl;
      std::cerr << "Error code: " << result << std::endl;
      return false;
    }
    pelem_ = pdoc_->RootElement();
    if (pelem_) {
    } else {
      std::cerr << "No root element found in XML" << std::endl;
    }
    return true;
  }

  /// Node name
  const char* GetName() const override {
    return pelem_->Name();
  }

  /// Attribute
  const char* GetAttribute(const char *name) const override {
    return pelem_->Attribute(name);
  }

  /// Walk the node
  void Walk(BasicNodeHandler& hnd) override {
    if (!pelem_)
      return;

    // Extract attributes
    const XMLAttribute* attr = pelem_->FirstAttribute();
    while (attr) {
      hnd.HandleAttribute(attr->Name(), attr->Value());
      attr = attr->Next();
    }

    // Process child elements
    XMLElement* child = pelem_->FirstChildElement();
    while (child) {
      TinyXML2Walker wlk(child);
      hnd.HandleSubnode(child->Name(), wlk);

      child = child->NextSiblingElement();
    }
  }

  /// Extract (pure) text value of the node
  std::string GetText() const override;
};


/// Make a default XML walker
std::unique_ptr<BasicTreeWalker>
MakeDefaultXMLWalker() {
  return std::make_unique<TinyXML2Walker>();
}

/// https://stackoverflow.com/a/43356508.
class TinyXML2TextWalker : public XMLVisitor
{
  std::string buf_;
public:
  /// Visit a text chunk
  bool Visit (const XMLText & txt) override
  {
    buf_ += txt.Value();
    return true;
  }
  /// Move out string
  std::string&& MoveOutString() { return std::move(buf_); }
};

std::string TinyXML2Walker::GetText() const {
  TinyXML2TextWalker prt;
  pelem_->Accept (&prt);
  return prt.MoveOutString();
}

}  // namespace mp

#if 0

// Claude draft

// Base callback class that users should inherit from
class MyCallback {
public:
  virtual ~MyCallback() = default;

  // Generic method called for each element with its name
  virtual void onElement(const std::string& elementName,
      const std::map<std::string, std::string>& attributes,
      const std::string& textContent) {
    std::cout << "Element: " << elementName << std::endl;
  }

  // Called when entering an element (before processing children)
  virtual void onElementStart(const std::string& elementName,
      const std::map<std::string, std::string>& attributes) {}

  // Called when exiting an element (after processing children)
  virtual void onElementEnd(const std::string& elementName) {}

  // Structure to hold mixed content information
  struct MixedContent {
    std::string plainText;
    std::vector<std::pair<std::string, std::string>> markupElements;
  };

  // Specific methods for each XML element type
  virtual void paramList(const std::map<std::string, std::string>& attributes) {}
  virtual void param(const std::map<std::string, std::string>& attributes) {}
  virtual void paramDescr(const std::string& text) {}
  virtual void paramTopic(const std::string& text) {}
  virtual void paramCategory(const std::string& text) {}
  virtual void paramValues(const std::map<std::string, std::string>& attributes) {}
  virtual void paramVal(const std::map<std::string, std::string>& attributes, const std::string& text) {}
  virtual void paramDefault(const std::string& text) {}
  virtual void paramNote(const std::string& text) {}
  virtual void paramNoteWithMarkup(const MixedContent& content) {
    // Default implementation just calls paramNote with plain text
    paramNote(content.plainText);
  }
};

// Main XML converter class
class MyXMLConverter {
private:
  MyCallback* callback;
  XMLDocument doc;

  // Get all text content from an element (including nested text and markup)
  std::string getTextContent(XMLElement* element) {
    std::string text;
    XMLNode* child = element->FirstChild();
    while (child) {
      if (child->ToText()) {
        const char* value = child->ToText()->Value();
        if (value) {
          text += value;
        }
      } else if (child->ToElement()) {
        // Handle inline markup elements (like <code>, <br/>)
        XMLElement* childElem = child->ToElement();
        std::string elemName = childElem->Name();

        if (elemName == "br") {
          text += "\n";
        } else if (elemName == "code") {
          const char* codeText = childElem->GetText();
          if (codeText) {
            text += codeText;
          }
        } else {
          // For other inline elements, just get their text
          const char* elemText = childElem->GetText();
          if (elemText) {
            text += elemText;
          }
        }
      }
      child = child->NextSibling();
    }
    // Trim leading/trailing whitespace but preserve internal formatting
    size_t start = text.find_first_not_of(" \t\n\r");
    size_t end = text.find_last_not_of(" \t\n\r");
    if (start != std::string::npos && end != std::string::npos) {
      return text.substr(start, end - start + 1);
    }
    return "";
  }

  // Get mixed content with markup information preserved
  MixedContent getMixedContent(XMLElement* element) {
    MixedContent content;
    std::string text;

    XMLNode* child = element->FirstChild();
    while (child) {
      if (child->ToText()) {
        const char* value = child->ToText()->Value();
        if (value) {
          text += value;
        }
      } else if (child->ToElement()) {
        XMLElement* childElem = child->ToElement();
        std::string elemName = childElem->Name();

        if (elemName == "br") {
          text += "\n";
          content.markupElements.push_back({"br", ""});
        } else if (elemName == "code") {
          const char* codeText = childElem->GetText();
          std::string codeContent = codeText ? codeText : "";
          text += codeContent;
          content.markupElements.push_back({"code", codeContent});
        } else {
          const char* elemText = childElem->GetText();
          std::string elemContent = elemText ? elemText : "";
          text += elemContent;
          content.markupElements.push_back({elemName, elemContent});
        }
      }
      child = child->NextSibling();
    }

    // Trim whitespace
    size_t start = text.find_first_not_of(" \t\n\r");
    size_t end = text.find_last_not_of(" \t\n\r");
    if (start != std::string::npos && end != std::string::npos) {
      content.plainText = text.substr(start, end - start + 1);
    } else {
      content.plainText = text;
    }

    return content;
  }

  // Recursive function to walk the XML tree
  void walkTree(XMLNode* node) {
    if (!node) return;

    XMLElement* element = node->ToElement();
    if (element) {
      std::string elementName = element->Name();

      // Extract attributes
      std::map<std::string, std::string> attributes;
      const XMLAttribute* attr = element->FirstAttribute();
      while (attr) {
        attributes[attr->Name()] = attr->Value();
        attr = attr->Next();
      }

      // Get text content
      std::string textContent = getTextContent(element);

      // Call element start callback
      callback->onElementStart(elementName, attributes);

      // Call specific method based on element name
      if (elementName == "paramList") {
        callback->paramList(attributes);
      } else if (elementName == "param") {
        callback->param(attributes);
      } else if (elementName == "paramDescr") {
        callback->paramDescr(textContent);
      } else if (elementName == "paramTopic") {
        callback->paramTopic(textContent);
      } else if (elementName == "paramCategory") {
        callback->paramCategory(textContent);
      } else if (elementName == "paramValues") {
        callback->paramValues(attributes);
      } else if (elementName == "paramVal") {
        callback->paramVal(attributes, textContent);
      } else if (elementName == "paramDefault") {
        callback->paramDefault(textContent);
      } else if (elementName == "paramNote") {
        // For paramNote, provide both plain text and mixed content
        MyCallback::MixedContent mixedContent = getMixedContent(element);
        callback->paramNoteWithMarkup(mixedContent);
        callback->paramNote(textContent);
      }

      // Call generic method
      callback->onElement(elementName, attributes, textContent);

      // Process child elements
      XMLNode* child = element->FirstChild();
      while (child) {
        if (child->ToElement()) {
          walkTree(child);
        }
        child = child->NextSibling();
      }

      // Call element end callback
      callback->onElementEnd(elementName);
    }
  }

public:
  MyXMLConverter(MyCallback* cb) : callback(cb) {}

  // Read and parse XML file
  bool readXMLFile(const std::string& filename) {
    XMLError result = doc.LoadFile(filename.c_str());
    if (result != XML_SUCCESS) {
      std::cerr << "Error loading XML file: " << filename << std::endl;
      std::cerr << "Error code: " << result << std::endl;
      return false;
    }
    return true;
  }

  // Walk the entire XML tree
  void convert() {
    XMLElement* root = doc.RootElement();
    if (root) {
      walkTree(root);
    } else {
      std::cerr << "No root element found in XML" << std::endl;
    }
  }
};

// Example usage with custom callback
class MyCustomCallback : public MyCallback {
private:
  int indentLevel = 0;

  void printIndent() {
    for (int i = 0; i < indentLevel; i++) {
      std::cout << "  ";
    }
  }

public:
  void onElementStart(const std::string& elementName,
      const std::map<std::string, std::string>& attributes) override {
    printIndent();
    std::cout << "-> " << elementName;
    if (!attributes.empty()) {
      std::cout << " [";
      bool first = true;
      for (const auto& attr : attributes) {
        if (!first) std::cout << ", ";
        std::cout << attr.first << "=\"" << attr.second << "\"";
        first = false;
      }
      std::cout << "]";
    }
    std::cout << std::endl;
    indentLevel++;
  }

  void onElementEnd(const std::string& elementName) override {
    indentLevel--;
  }

  void param(const std::map<std::string, std::string>& attributes) override {
    printIndent();
    std::cout << "*** Processing parameter: " << attributes.at("name") << std::endl;
  }

  void paramDescr(const std::string& text) override {
    printIndent();
    std::cout << "Description: " << text << std::endl;
  }

  void paramVal(const std::map<std::string, std::string>& attributes,
      const std::string& text) override {
    printIndent();
    std::cout << "Value " << attributes.at("value") << ": " << text << std::endl;
  }

  void paramNote(const std::string& text) override {
    printIndent();
    std::cout << "Note: " << text.substr(0, 100) << "..." << std::endl;
  }

  void paramNoteWithMarkup(const MyCallback::MixedContent& content) override {
    printIndent();
    std::cout << "Note with markup (plain text): " << std::endl;
    printIndent();
    std::cout << content.plainText.substr(0, 200) << "..." << std::endl;

    if (!content.markupElements.empty()) {
      printIndent();
      std::cout << "Contains " << content.markupElements.size() << " markup elements:" << std::endl;
      for (const auto& elem : content.markupElements) {
        printIndent();
        std::cout << "  <" << elem.first << ">";
        if (!elem.second.empty()) {
          std::cout << elem.second;
        }
        std::cout << std::endl;
      }
    }
  }
};

// Main function demonstrating usage
int main() {
  MyCustomCallback callback;
  MyXMLConverter converter(&callback);

  if (converter.readXMLFile("params.xml")) {
    std::cout << "=== Converting XML ===" << std::endl;
    converter.convert();
    std::cout << "=== Conversion complete ===" << std::endl;
  } else {
    std::cerr << "Failed to read XML file" << std::endl;
    return 1;
  }

  return 0;
}


#endif // ...
