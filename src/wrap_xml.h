#ifndef WRAP_XML_H
#define WRAP_XML_H

#include <string>
#include <memory>

namespace mp {

class BasicNodeHandler;

class BasicTreeWalker;


/// This is specialized by the user
/// to handle a specific node
class BasicNodeHandler {
public:
  /// Destroy
  virtual ~BasicNodeHandler() { }
  /// Handle atribute
  virtual void HandleAttribute(
      const char* name, const char* value) { }
  /// Handle subnode (element)
  virtual void HandleSubnode(
      const char* name, BasicTreeWalker& walker) { }
};


/// This stores either the whole tree,
/// or a subnode.
class BasicTreeWalker {
public:
  /// Destroy
  virtual ~BasicTreeWalker() { }
  /// Read file (for the root only)
  virtual bool ReadFile(const char* name) = 0;
  /// Node name
  virtual const char* GetName() const = 0;
  /// Attribute with given name, if any
  virtual const char* GetAttribute(const char* name) const = 0;
  /// Walk the node's children
  virtual void Walk(BasicNodeHandler& hnd) = 0;
  /// Extract (unformatted) text value of the node
  virtual std::string GetText() const = 0;
};


/// Make a default XML walker
std::unique_ptr<BasicTreeWalker>
MakeDefaultXMLWalker();

}  // namespace mp

#endif // WRAP_XML_H
