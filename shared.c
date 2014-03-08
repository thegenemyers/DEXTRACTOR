typedef unsigned long long uint64;
typedef unsigned int       uint32;
typedef unsigned short     uint16;
typedef unsigned char      uint8;

static char *Program_Name;

static void *Guarded_Alloc(void *p, size_t size)
{ p = realloc(p,size);
  if (p == NULL)
    { fprintf(stderr,"%s: Out of memory\n",Program_Name);
      exit (1);
    }
  return (p);
}

static FILE *Guarded_Fopen(char *name, char *io)
{ FILE *f = fopen(name,io);
  if (f == NULL)
    { fprintf(stderr,"%s: Cannot open file %s\n",Program_Name,name);
      exit (1);
    }
  return (f);
}

#ifdef DEXTRACT

  static char *Guarded_Strdup(char *name)
  { char *s = strdup(name);
    if (s == NULL)
      { fprintf(stderr,"%s: Out of memory\n",Program_Name);
        exit (1);
      }
    return (s);
  }

#else

#define COMPRESSED_LEN(len)  (((len)+3) >> 2)

#if defined(DEXTA) || defined(DEXQV)

  //  Compress read s[0..len-1] into 2 bit format (in place)

  void Compress_Read(int len, char *s)

  { static char number[128] =
      { 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 3, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 3, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
      };

    int    i, c, d;
    uint8 *s0, *s1, *s2, *s3;

    s0 = (uint8 *) s;
    s1 = s0+1;
    s2 = s1+1;
    s3 = s2+1;

    c = s1[len];
    d = s2[len];
    s0[len] = s1[len] = s2[len] = 0;

    for (i = 0; i < len; i += 4)
      *s++ = (number[s0[i]] << 6) | (number[s1[i]] << 4) | (number[s2[i]] << 2) | number[s3[i]];

    s1[len] = c;
    s2[len] = d;
  }

#else  // defined(UNDEXTA) || defined(UNDEXQV)

  //  Uncompress 2-bit packed string in s, in place

  void Uncompress_Read(int len, char *s)

  { static char letter[4] = { 'A', 'C', 'G', 'T' };

    int   i, tlen, byte;
    char *s0, *s1, *s2, *s3;
    char *t;

    s0 = s;
    s1 = s0+1;
    s2 = s1+1;
    s3 = s2+1;

    tlen = (len-1)/4;

    t = s+tlen;
    for (i = tlen*4; i >= 0; i -= 4)
      { byte = *t--;
        s0[i] = letter[((byte >> 6) & 0x3)];
        s1[i] = letter[((byte >> 4) & 0x3)];
        s2[i] = letter[((byte >> 2) & 0x3)];
        s3[i] = letter[(byte & 0x3)];
      }
    s[len] = '\0';
  }

  static void flip_short(void *w)
  { uint8 *v = (uint8 *) w;
    uint8  x;

    x    = v[0];
    v[0] = v[1];
    v[1] = x;
  }

  static void flip_long(void *w)
  { uint8 *v = (uint8 *) w;
    uint8  x;

    x    = v[0];
    v[0] = v[3];
    v[3] = x;
    x    = v[1];
    v[1] = v[2];
    v[2] = x;
  }

#endif
#endif
