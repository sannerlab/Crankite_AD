// structure for full atoms rotameric side chains
struct _AASCRot {
  int nbRot;
  int nbAtoms;
  int *atypes;      /* list of nbAtoms integers describing atom types */
  char *name;       /* residue name */
  double *charges;  /* list of nbAtoms doubles holding atomic charges */
  char *atnames;
  double **coords;  /* [nbRot][natoms x 3] double holding atomic coordinates for rotamers */
};

struct _ARG {
  int nbRot;
  int nbAtoms;
  int atypes[11];
  double charges[11];
  char atnames[11*5];
  double coords[81][11][3];
};

struct _ASN {
  int nbRot;
  int nbAtoms;
  int atypes[5];
  double charges[5];
  char atnames[5*5];
  double coords[18][5][3];
};

struct _ASP {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  char atnames[3*5];
  double coords[9][3][3];
};

struct _CYS {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  char atnames[2*5];
  double coords[3][2][3];
};

struct _GLN {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[36][6][3];
};

struct _GLU {
  int nbRot;
  int nbAtoms;
  int atypes[4];
  double charges[4];
  char atnames[4*5];
  double coords[27][4][3];
};

struct _HIS {
  int nbRot;
  int nbAtoms;
  int atypes[5];
  double charges[5];
  char atnames[5*5];
  double coords[9][5][3];
};

struct _HID {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[9][6][3];
};

struct _HIE {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[9][6][3];
};

struct _ILE {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  char atnames[3*5];
  double coords[9][3][3];
};

struct _LEU {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  char atnames[3*5];
  double coords[9][3][3];
};

struct _LYS {
  int nbRot;
  int nbAtoms;
  int atypes[7];
  double charges[7];
  char atnames[7*5];
  double coords[81][7][3];
};

struct _MET {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  char atnames[3*5];
  double coords[27][3][3];
};

struct _PHE {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[6][6][3];
};

struct _PRO {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  char atnames[2*5];
  double coords[2][2][3];
};

struct _SER {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  char atnames[2*5];
  double coords[3][2][3];
};

struct _THR {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  char atnames[3*5];
  double coords[3][3][3];
};

struct _TRP {
  int nbRot;
  int nbAtoms;
  int atypes[10];
  double charges[10]; 
  char atnames[10*5];
  double coords[9][10][3];
};

struct _TYR {
  int nbRot;
  int nbAtoms;
  int atypes[8];
  double charges[8];
  char atnames[8*5];
  double coords[6][8][3];
};

struct _VAL {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  char atnames[2*5];
  double coords[3][2][3];
};

struct _Af1 {
  int nbRot;
  int nbAtoms;
  int atypes[1];
  double charges[1];
  char atnames[1*5];
  double coords[3][1][3];
};

struct _Af2 {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  char atnames[2*5];
  double coords[3][2][3];
};

struct _Af3 {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  char atnames[3*5];
  double coords[1][3][3];
};

struct _Bf1 {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  char atnames[2*5];
  double coords[9][2][3];
};

struct _Bf2 {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  char atnames[3*5];
  double coords[9][3][3];
};

struct _Bf3 {
  int nbRot;
  int nbAtoms;
  int atypes[4];
  double charges[4];
  char atnames[4*5];
  double coords[3][4][3];
};

struct _Cf3 {
  int nbRot;
  int nbAtoms;
  int atypes[5];
  double charges[5];
  char atnames[5*5];
  double coords[18][5][3];
};

struct _Ff3 {
  int nbRot;
  int nbAtoms;
  int atypes[9];
  double charges[9];
  char atnames[9*5];
  double coords[2][9][3];
};

struct _If2 {
  int nbRot;
  int nbAtoms;
  int atypes[4];
  double charges[4];
  char atnames[4*5];
  double coords[9][4][3];
};

struct _If3b {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[6][6][3];
};

struct _If3g {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[6][6][3];
};

struct _Lf2 {
  int nbRot;
  int nbAtoms;
  int atypes[4];
  double charges[4];
  char atnames[4*5];
  double coords[6][4][3];
};

struct _Lf3r {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[4][6][3];
};

struct _Lf3s {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[4][6][3];
};

struct _Lf6 {
  int nbRot;
  int nbAtoms;
  int atypes[9];
  double charges[9];
  char atnames[9*5];
  double coords[4][9][3];
};

struct _Mf3 {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[18][6][3];
};

struct _Vf1 {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  char atnames[3*5];
  double coords[2][3][3];
};

struct _Vf3r {
  int nbRot;
  int nbAtoms;
  int atypes[5];
  double charges[5];
  char atnames[3*5];
  double coords[3][5][3];
};

struct _Vf3s {
  int nbRot;
  int nbAtoms;
  int atypes[5];
  double charges[5];
  char atnames[5*5];
  double coords[2][5][3];
};

struct _Vf6 {
  int nbRot;
  int nbAtoms;
  int atypes[8];
  double charges[8];
  char atnames[6*5];
  double coords[2][8][3];
};

struct _SEP {
  int nbRot;
  int nbAtoms;
  int atypes[5];
  double charges[5];
  char atnames[5*5];
  double coords[9][5][3];
};

struct _TPO {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  char atnames[6*5];
  double coords[27][6][3];
};

struct _PTR {
  int nbRot;
  int nbAtoms;
  int atypes[11];
  double charges[11];
  char atnames[11*5];
  double coords[81][11][3];
};



extern int getSideChainTemplateIndexFromName(char *str);
extern int getSideChainTemplateIndexFromIDchar(char id);
extern void intialize_AASCRotTable();
extern struct _AASCRot *_AASCRotTable;
