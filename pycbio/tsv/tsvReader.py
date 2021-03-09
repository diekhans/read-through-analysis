# Copyright 2006-2012 Mark Diekhans
"""TSV reading classes"""
import sys
import csv
from pycbio.sys import fileOps
from pycbio.tsv.tsvRow import TsvRow
from pycbio.tsv import TsvError

csv.field_size_limit(sys.maxsize)

# FIXME:  pass owndership of row to Row instead of having Row inherit from list
# FIXME: make rowClass (rename rowFactory) so that it can be able to construct any
#        object.  Move parsing of columns outside of tsvRow (good for peewee)
# FIXME: create error with file name/line number
# FIMXE: carefully consider naming of functions in Row, as they could
#  conflict with fields.  Maybe put in a base-class??
# FIXME: put in same module as TSV
# FIXME: TsvError and preserving the traceback is a pain
# FIXME: is colMap needed any more???
# FIXME: need to add write stuff. (see GeneCheck), where str() is called on
#        alfor all column ty
# FIXME: add type mapping functions that gets column name
# FIXME: rename  typeMap -> colTypes
# FIXME: make a column object.
# FIXME: document how colName mapping and type mappings works together
# FIXME: check if column is a valid python field name
# FIXME: switch to row derived from namedtuple
# FIXME: rowClass interface is hacky.  It could be a keyword/value and not have to do column lookup.
# FIXME: maybe build on csv.Reader class and keep less of our own crap (although column stuff is nice)
#        however, csv.reader has got us into trouble because of print formatted files with quotes in data
# FIXME: add default for column not in file to typemap

# typeMap converter for str types were empty represents None
strOrNoneType = (lambda v: None if (v == "") else v,
                 lambda v: "" if (v is None) else v)

# typeMap converter for int types were empty represents None
intOrNoneType = (lambda v: None if (v == "") else int(v),
                 lambda v: "" if (v is None) else str(v))

floatOrNoneType = (lambda v: None if (v == "") else float(v),
                   lambda v: "" if (v is None) else str(v))


class printf_basic_dialect(csv.Dialect):
    """Describes the usual properties for TSV files generated by printf, etc.  Common
    in bioinformatics.  Quotes can be included in data, etc.
    """
    delimiter = '\t'
    quotechar = ''
    doublequote = False
    skipinitialspace = False
    lineterminator = '\n'
    quoting = csv.QUOTE_NONE


class TsvReader(object):
    """Class for reading TSV files.  Reads header and builds column name to
    column index map.  After a next, object contains a row and each column
    becomes a field name.  It is also indexable by column name or int index.
    Columns can be automatically type converted by column name.  This can also
    read from a dbapi cursor object (must set allowEmpty to true)

    If the first character of the header is '#', it is skipped.
    """

    def _readRow(self):
        "read the next row, returning None on EOF"
        if self.reader is None:
            return None
        try:
            row = next(self.reader)
        except Exception as ex:
            self.close()
            if isinstance(ex, StopIteration):
                return None
            else:
                raise
        self.lineNum = self.reader.line_num
        return row

    def _readHeader(self, allowEmpty):
        row = self._readRow()
        if row is None:
            if not allowEmpty:
                raise TsvError("empty TSV file", reader=self)
        else:
            if (len(row) > 0) and row[0].startswith('#'):
                row[0] = row[0][1:]
            self._setupColumns(row)

    def _setupColumns(self, columns):
        # n.b. columns could be passed in from client, must copy
        i = 0
        for col in columns:
            if self.columnNameMapper is not None:
                self.extColumns.append(col)
                col = self.columnNameMapper(col)
            self.columns.append(col)
            if col in self.colMap:
                raise TsvError("Duplicate column name: {}".format(col))
            self.colMap[col] = i
            i += 1

    def _initColTypes(self, typeMap, defaultColType):
        "save col types as column indexed list"
        if typeMap is not None:
            # build from type map
            self.colTypes = []
            for col in self.columns:
                self.colTypes.append(typeMap.get(col, defaultColType))
        elif defaultColType is not None:
            # fill in colTypes from default
            self.colTypes = []
            for i in range(len(self.columns)):
                self.colTypes.append(defaultColType)

    def __init__(self, fileName, rowClass=None, typeMap=None, defaultColType=None, columns=None, columnNameMapper=None,
                 ignoreExtraCols=False, inFh=None, allowEmpty=False, dialect=csv.excel_tab,
                 encoding=None, errors=None):
        """Open TSV file and read header into object.  Removes leading # from
        UCSC header.

        fileName - name of file, opened unless inFh is specified
        rowClass - class or class factory function to use for a row. Must take
            TsvReader and list of string values of columns.
        typeMap - if specified, it maps column names to the type objects to
            use to convert the column.  Unspecified columns will not be
            converted. Key is the column name, value can be either a type
            or a tuple of (parseFunc, formatFunc).  If a type is use,
            str() is used to convert to a printable value.
        defaultColType - if specified, type of unspecified columns
        columns - if specified, the column names to use.  The header
            should not be in the file.
        columnNameMapper - function to map column names to the internal name.
        ignoreExtraCols - should extra columns be ignored?
        inFh - If not None, this is used as the open file, rather than
          opening it.  Closed when the end of file is reached.
        allowEmpty - an empty input results in an EOF rather than an error.
        dialect - a csv dialect object or name.
        """
        self.columns = []
        # external column name, before mapping with columnNameMapper
        if columnNameMapper is None:
            self.extColumns = self.columns
        else:
            self.extColumns = []
        self.colMap = {}
        self.fileName = fileName
        self.lineNum = 0
        self.rowClass = rowClass
        if rowClass is None:
            self.rowClass = TsvRow
        self.columnNameMapper = columnNameMapper
        self.colTypes = None
        self.ignoreExtraCols = ignoreExtraCols
        if inFh is not None:
            self.inFh = inFh
        else:
            self.inFh = fileOps.opengz(fileName, encoding=encoding, errors=errors)
        try:
            self.reader = csv.reader(self.inFh, dialect=dialect)
            if columns:
                self._setupColumns(columns)
            else:
                self._readHeader(allowEmpty)
            self._initColTypes(typeMap, defaultColType)
        except Exception:
            self.close()
            raise

    def close(self):
        """close file, keeping column map around.  Called automatically
        by iter on EOF."""
        if self.inFh is not None:
            self.inFh.close()
            self.inFh = None
            self.reader = None

    def __iter__(self):
        return self

    def __next__(self):
        try:
            row = self._readRow()
        except Exception as ex:
            raise TsvError("Error reading TSV row", self) from ex
        if row is None:
            raise StopIteration
        if ((self.ignoreExtraCols and (len(row) < len(self.columns))) or ((not self.ignoreExtraCols) and (len(row) != len(self.columns)))):
            # FIXME: will hang: self.close()
            raise TsvError("row has {} columns, expected {}".format(len(row), len(self.columns)), reader=self)
        try:
            return self.rowClass(self, row)
        except Exception as ex:
            raise TsvError("Error converting TSV row to object", self) from ex
