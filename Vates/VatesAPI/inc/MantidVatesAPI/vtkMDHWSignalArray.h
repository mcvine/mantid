/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMDHWSignalArray.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkMDHWSignalArray - Map native Exodus II results arrays
// into the vtkDataArray interface.
//
// .SECTION Description
// Map native Exodus II results arrays into the vtkDataArray interface. Use
// the vtkCPExodusIIInSituReader to read an Exodus II file's data into this
// structure.

#ifndef vtkMDHWSignalArray_h
#define vtkMDHWSignalArray_h

#include "vtkMappedDataArray.h"

#include "vtkTypeTemplate.h"  // For templated vtkObject API
#include "vtkObjectFactory.h" // for vtkStandardNewMacro
#include "vtkIdList.h"
#include "vtkVariant.h"
#include "vtkVariantCast.h"

#include "MantidDataObjects/MDHistoWorkspaceIterator.h"

namespace Mantid {
namespace VATES {

template <class Scalar>
class vtkMDHWSignalArray : public vtkTypeTemplate<vtkMDHWSignalArray<Scalar>,
                                                  vtkMappedDataArray<Scalar>> {
public:
  vtkMappedDataArrayNewInstanceMacro(
      vtkMDHWSignalArray<Scalar>) static vtkMDHWSignalArray *New();
  void PrintSelf(ostream &os, vtkIndent indent) override;

  void InitializeArray(
      std::unique_ptr<Mantid::DataObjects::MDHistoWorkspaceIterator> iterator,
      std::size_t offset, vtkIdType size);

  // Reimplemented virtuals -- see superclasses for descriptions:
  void Initialize() override;
  void GetTuples(vtkIdList *ptIds, vtkAbstractArray *output) override;
  void GetTuples(vtkIdType p1, vtkIdType p2, vtkAbstractArray *output) override;
  void Squeeze() override;
  vtkArrayIterator *NewIterator() override;
  vtkIdType LookupValue(vtkVariant value) override;
  void LookupValue(vtkVariant value, vtkIdList *ids) override;
  vtkVariant GetVariantValue(vtkIdType idx) override;
  void ClearLookup() override;
  double *GetTuple(vtkIdType i) override;
  void GetTuple(vtkIdType i, double *tuple) override;
  vtkIdType LookupTypedValue(Scalar value) override;
  void LookupTypedValue(Scalar value, vtkIdList *ids) override;
  Scalar GetValue(vtkIdType idx) override;
  Scalar &GetValueReference(vtkIdType idx) override;
  void GetTupleValue(vtkIdType idx, Scalar *t) override;

  // Description:
  // This container is read only -- this method does nothing but print a
  // warning.
  int Allocate(vtkIdType sz, vtkIdType ext) override;
  int Resize(vtkIdType numTuples) override;
  void SetNumberOfTuples(vtkIdType number) override;
  void SetTuple(vtkIdType i, vtkIdType j, vtkAbstractArray *source) override;
  void SetTuple(vtkIdType i, const float *source) override;
  void SetTuple(vtkIdType i, const double *source) override;
  void InsertTuple(vtkIdType i, vtkIdType j, vtkAbstractArray *source) override;
  void InsertTuple(vtkIdType i, const float *source) override;
  void InsertTuple(vtkIdType i, const double *source) override;
  void InsertTuples(vtkIdList *dstIds, vtkIdList *srcIds,
                    vtkAbstractArray *source) override;
  void InsertTuples(vtkIdType dstStart, vtkIdType n, vtkIdType srcStart,
                    vtkAbstractArray *source) override;
  vtkIdType InsertNextTuple(vtkIdType j, vtkAbstractArray *source) override;
  vtkIdType InsertNextTuple(const float *source) override;
  vtkIdType InsertNextTuple(const double *source) override;
  void DeepCopy(vtkAbstractArray *aa) override;
  void DeepCopy(vtkDataArray *da) override;
  void InterpolateTuple(vtkIdType i, vtkIdList *ptIndices,
                        vtkAbstractArray *source, double *weights) override;
  void InterpolateTuple(vtkIdType i, vtkIdType id1, vtkAbstractArray *source1,
                        vtkIdType id2, vtkAbstractArray *source2,
                        double t) override;
  void SetVariantValue(vtkIdType idx, vtkVariant value) override;
  void RemoveTuple(vtkIdType id) override;
  void RemoveFirstTuple() override;
  void RemoveLastTuple() override;
  void SetTupleValue(vtkIdType i, const Scalar *t) override;
  void InsertTupleValue(vtkIdType i, const Scalar *t) override;
  vtkIdType InsertNextTupleValue(const Scalar *t) override;
  void SetValue(vtkIdType idx, Scalar value) override;
  vtkIdType InsertNextValue(Scalar v) override;
  void InsertVariantValue(vtkIdType idx, vtkVariant value) override;
  void InsertValue(vtkIdType idx, Scalar v) override;

protected:
  vtkMDHWSignalArray();
  ~vtkMDHWSignalArray() override;

private:
  vtkIdType Lookup(const Scalar &val, vtkIdType startIndex);
  vtkMDHWSignalArray(const vtkMDHWSignalArray &); // Not implemented.
  void operator=(const vtkMDHWSignalArray &);     // Not implemented.

  std::size_t m_offset;
  std::unique_ptr<Mantid::DataObjects::MDHistoWorkspaceIterator> m_iterator;
  Scalar m_temporaryTuple[1];
};

//! @cond Doxygen_Suppress
//------------------------------------------------------------------------------
/// Can't use vtkStandardNewMacro on a templated class.
/// create and return a pointer to a vtkMDHWSignalArray<Scalar> object.
template <class Scalar>
vtkMDHWSignalArray<Scalar> *vtkMDHWSignalArray<Scalar>::New() {
  VTK_STANDARD_NEW_BODY(vtkMDHWSignalArray<Scalar>)
}
//! @endcond
//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::PrintSelf(ostream &os, vtkIndent indent) {
  this->vtkMDHWSignalArray<Scalar>::Superclass::PrintSelf(os, indent);
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InitializeArray(
    std::unique_ptr<Mantid::DataObjects::MDHistoWorkspaceIterator> iterator,
    std::size_t offset, vtkIdType size) {
  this->MaxId = size - 1;
  this->Size = size;
  this->NumberOfComponents = 1;
  this->m_iterator = std::move(iterator);
  this->m_offset = offset;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkMDHWSignalArray<Scalar>::Initialize() {
  this->MaxId = -1;
  this->Size = 0;
  this->NumberOfComponents = 1;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::GetTuples(vtkIdList *ptIds,
                                           vtkAbstractArray *output) {
  vtkDataArray *da = vtkDataArray::FastDownCast(output);
  if (!da) {
    vtkWarningMacro(<< "Input is not a vtkDataArray");
    return;
  }

  if (da->GetNumberOfComponents() != this->GetNumberOfComponents()) {
    vtkWarningMacro(<< "Incorrect number of components in input array.");
    return;
  }

  const vtkIdType numPoints = ptIds->GetNumberOfIds();
  for (vtkIdType i = 0; i < numPoints; ++i) {
    da->SetTuple(i, this->GetTuple(ptIds->GetId(i)));
  }
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::GetTuples(vtkIdType p1, vtkIdType p2,
                                           vtkAbstractArray *output) {
  vtkDataArray *da = vtkDataArray::FastDownCast(output);
  if (!da) {
    vtkErrorMacro(<< "Input is not a vtkDataArray");
    return;
  }

  if (da->GetNumberOfComponents() != this->GetNumberOfComponents()) {
    vtkErrorMacro(<< "Incorrect number of components in input array.");
    return;
  }

  for (vtkIdType daTupleId = 0; p1 <= p2; ++p1) {
    da->SetTuple(daTupleId++, this->GetTuple(p1));
  }
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkMDHWSignalArray<Scalar>::Squeeze() {
  // noop
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkArrayIterator *vtkMDHWSignalArray<Scalar>::NewIterator() {
  vtkErrorMacro(<< "Not implemented.");
  return NULL;
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkIdType vtkMDHWSignalArray<Scalar>::LookupValue(vtkVariant value) {
  bool valid = true;
  Scalar val = vtkVariantCast<Scalar>(value, &valid);
  if (valid) {
    return this->Lookup(val, 0);
  }
  return -1;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::LookupValue(vtkVariant value, vtkIdList *ids) {
  bool valid = true;
  Scalar val = vtkVariantCast<Scalar>(value, &valid);
  ids->Reset();
  if (valid)
    {
    vtkIdType index = 0;
    while ((index = this->Lookup(val, index)) >= 0)
      {
      ids->InsertNextId(index);
      ++index;
      }
    }
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkVariant vtkMDHWSignalArray<Scalar>::GetVariantValue(vtkIdType idx) {
  return vtkVariant(this->GetValue(idx));
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkMDHWSignalArray<Scalar>::ClearLookup() {
  // no-op, no fast lookup implemented.
}

//------------------------------------------------------------------------------
template <class Scalar>
double *vtkMDHWSignalArray<Scalar>::GetTuple(vtkIdType i) {
  m_iterator->jumpTo(m_offset + i);
  m_temporaryTuple[0] = m_iterator->getNormalizedSignal();
  return &m_temporaryTuple[0];
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::GetTuple(vtkIdType i, double *tuple) {
  m_iterator->jumpTo(m_offset + i);
  tuple[0] = m_iterator->getNormalizedSignal();
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkIdType vtkMDHWSignalArray<Scalar>::LookupTypedValue(Scalar value) {
  return this->Lookup(value, 0);
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::LookupTypedValue(Scalar value,
                                                  vtkIdList *ids) {
  ids->Reset();
  vtkIdType index = 0;
  while ((index = this->Lookup(value, index)) >= 0)
    {
    ids->InsertNextId(index);
    ++index;
    }
}

template <class Scalar>
vtkIdType vtkMDHWSignalArray<Scalar>::Lookup(const Scalar &val,
                                             vtkIdType index) {
  while (index <= this->MaxId) {
    if (this->GetValueReference(index) == val) {
      return index;
    }
    index++;
  }
  return -1;
}

//------------------------------------------------------------------------------
template <class Scalar>
Scalar vtkMDHWSignalArray<Scalar>::GetValue(vtkIdType idx) {
  m_iterator->jumpTo(m_offset + idx);
  return m_iterator->getNormalizedSignal();
}
//------------------------------------------------------------------------------
template <class Scalar>
Scalar &vtkMDHWSignalArray<Scalar>::GetValueReference(vtkIdType idx) {
  m_iterator->jumpTo(m_offset + idx);
  m_temporaryTuple[0] = m_iterator->getNormalizedSignal();
  return m_temporaryTuple[0];
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::GetTupleValue(vtkIdType tupleId,
                                               Scalar *tuple) {
  m_iterator->jumpTo(m_offset + tupleId);
  tuple[0] = m_iterator->getNormalizedSignal();
}

//------------------------------------------------------------------------------
template <class Scalar>
int vtkMDHWSignalArray<Scalar>::Allocate(vtkIdType, vtkIdType) {
  vtkErrorMacro("Read only container.") return 0;
}

//------------------------------------------------------------------------------
template <class Scalar> int vtkMDHWSignalArray<Scalar>::Resize(vtkIdType) {
  vtkErrorMacro("Read only container.") return 0;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::SetNumberOfTuples(vtkIdType) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::SetTuple(vtkIdType, vtkIdType,
                                          vtkAbstractArray *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::SetTuple(vtkIdType, const float *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::SetTuple(vtkIdType, const double *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InsertTuple(vtkIdType, vtkIdType,
                                             vtkAbstractArray *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InsertTuple(vtkIdType, const float *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InsertTuple(vtkIdType, const double *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InsertTuples(vtkIdList *, vtkIdList *,
                                              vtkAbstractArray *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InsertTuples(vtkIdType, vtkIdType, vtkIdType,
                                              vtkAbstractArray *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkIdType vtkMDHWSignalArray<Scalar>::InsertNextTuple(vtkIdType,
                                                      vtkAbstractArray *) {
  vtkErrorMacro("Read only container.") return -1;
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkIdType vtkMDHWSignalArray<Scalar>::InsertNextTuple(const float *) {

  vtkErrorMacro("Read only container.") return -1;
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkIdType vtkMDHWSignalArray<Scalar>::InsertNextTuple(const double *) {
  vtkErrorMacro("Read only container.") return -1;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::DeepCopy(vtkAbstractArray *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::DeepCopy(vtkDataArray *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InterpolateTuple(vtkIdType, vtkIdList *,
                                                  vtkAbstractArray *,
                                                  double *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InterpolateTuple(vtkIdType, vtkIdType,
                                                  vtkAbstractArray *, vtkIdType,
                                                  vtkAbstractArray *, double) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::SetVariantValue(vtkIdType, vtkVariant) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::RemoveTuple(vtkIdType) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkMDHWSignalArray<Scalar>::RemoveFirstTuple() {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkMDHWSignalArray<Scalar>::RemoveLastTuple() {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::SetTupleValue(vtkIdType, const Scalar *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InsertTupleValue(vtkIdType, const Scalar *) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkIdType vtkMDHWSignalArray<Scalar>::InsertNextTupleValue(const Scalar *) {
  vtkErrorMacro("Read only container.") return -1;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::SetValue(vtkIdType, Scalar) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkIdType vtkMDHWSignalArray<Scalar>::InsertNextValue(Scalar) {
  vtkErrorMacro("Read only container.") return -1;
}

//------------------------------------------------------------------------------
template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InsertValue(vtkIdType, Scalar) {
  vtkErrorMacro("Read only container.") return;
}

template <class Scalar>
void vtkMDHWSignalArray<Scalar>::InsertVariantValue(vtkIdType, vtkVariant) {
  vtkErrorMacro("Read only container.") return;
}

//------------------------------------------------------------------------------
template <class Scalar>
vtkMDHWSignalArray<Scalar>::vtkMDHWSignalArray()
    : m_offset(0) {}

//------------------------------------------------------------------------------
template <class Scalar> vtkMDHWSignalArray<Scalar>::~vtkMDHWSignalArray() {}
}
}

#endif // vtkMDHWSignalArray_h

// VTK-HeaderTest-Exclude: vtkMDHWSignalArray.h
