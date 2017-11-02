// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#include "ReactionProcessor.h"
#include "Node.h"

namespace nla3d {

ReactionProcessor::ReactionProcessor(FEStorage *st) : PostProcessor(st) {
  name ="ReactionProcessor";
}

ReactionProcessor::ReactionProcessor(FEStorage *st, std::string _filename) : PostProcessor(st) {
  name ="ReactionProcessor";
  filename = _filename;
}

void ReactionProcessor::pre() {
  if (nodes.size() == 0) {
    LOG(WARNING) << "Can't work. No nodes. Processor name = " << name;
    return;
  }

  // define Dofs by the following rule:
  // 1. If dofs is not empty than user alreade define which dof use to calculate reaction
  // 2. If dofs.size() == 0 , select dof thich are constrained in current node set.
  // if different nodes has different constrained dofs - error
  if (dofs.size() == 0) {
    for (uint16 d = 0; d < Dof::numberOfDofTypes; d++) {
      Dof::dofType t = static_cast<Dof::dofType> (d);
      if (storage->isNodeDofUsed(nodes[0], t) && storage->getNodeDof(nodes[0], t)->isConstrained) {
        dofs.push_back(t);
      }
    }
    if (dofs.size() == 0) {
      LOG(WARNING) << "Can't select dofs for reactions (they are unconstrained)";
    }
    // check that others nodes have the same constrained dofs
    bool sameDofsConstrained = true;
    uint32 n = 0;
    while (sameDofsConstrained && n < nodes.size()) {
      for (uint16 d = 0; d < dofs.size(); d++) {
        if (!storage->isNodeDofUsed(nodes[n], dofs[d]) && !storage->getNodeDof(nodes[n], dofs[d])->isConstrained) {
          LOG(FATAL) << "Different dofs are constrained in the node set. Autochoosing of dofs is failed";
          sameDofsConstrained = false;
          break;
        }
      }
      n++;
    }
  } else {
    // in case of user defined dofs just check that this dofs are used in the FE
    bool dofsUsed = true;
    uint32 n = 0;
    while (dofsUsed && n < nodes.size()) {
      for (uint16 d = 0; d < dofs.size(); d++) {
        if (!storage->isNodeDofUsed(nodes[n], dofs[d])) {
          LOG(ERROR) << "some dofs are not used in FE calculations";
          dofsUsed = false;
          break;
        }
      }
      n++;
    }
  }
  if (filename.length() > 0) {
    std::ofstream file(filename.c_str(),std::ios::trunc);
    if (!file) {
      LOG(WARNING) << "Can't create a file with name " << filename;
      return;
    }
    for (uint16 d = 0; d < dofs.size(); d++) {
      file << "\t" << Dof::dofTypeLabels[dofs[d]];
    }
    file << std::endl;
    file << "0";
    for (uint16 d = 0; d < dofs.size(); d++) {
      file << "\t" << 0.0;
    }
    file << std::endl;
    file.close();
  }
  sumOfDofsReactions.assign(dofs.size(), std::vector<double> ());
  for (uint16 d = 0; d < dofs.size(); d++) {
    //sumOfDofsReactions[d].reserve(qLoadstep+1);
    sumOfDofsReactions[d].push_back(0.0);
  }
}

void ReactionProcessor::process (uint16 curLoadstep) {
  std::vector<double> reactions;
  reactions.assign(dofs.size(), 0.0);

  for (uint32 n = 0; n < nodes.size(); n++) {
    for (uint16 d = 0; d < dofs.size(); d++) {
      reactions[d] += storage->getReaction(nodes[n],dofs[d]);
    }
  }

  if (filename.length() > 0) {
    std::ofstream file(filename.c_str(),std::ios::app);
    if (!file) {
      LOG(WARNING) << "Can't open a file with name " << filename;
      return;
    }
    file << curLoadstep;
    for (uint16 d = 0; d < dofs.size(); d++) {
      file << "\t" << reactions[d];
    }
    file << std::endl;
    file.close();
  }
  for (uint16 d = 0; d < dofs.size(); d++) {
    sumOfDofsReactions[d].push_back(reactions[d]);
  }
}


void ReactionProcessor::post (uint16 curLoadstep) {
}


std::vector<double> ReactionProcessor::getReactions (Dof::dofType dof) {
  for (uint16 d = 0; d < dofs.size(); d++) {
    if (dofs[d] == dof) {
      return sumOfDofsReactions[d];
    }
  }
  LOG(ERROR) << "There is not results for DoF = " << Dof::dofTypeLabels[dof];
  return std::vector<double> ();
}

} // namespace nla3d
