// vim: sw=2 ts=2 expandtab smartindent

#define WASM_EXPORT __attribute__((visibility("default")))
#define ARR_LEN(arr) (sizeof(arr) / sizeof((arr)[0]))

/* get a grip, man!
 * ooh, gettin handsy, hm?
 * get a hold of yourself!
 * gotta hand it to ya ...
 * need a hand with that?
 * so fun, it's disarming!
 */

#define MATH_PI (3.141592653589793)
#define MATH_TAU (6.283185307179586)
#define TICK_SECOND 60
#include <stdint.h>
#include "font.h"

/* -- BEGIN WASM -- */
extern void print(float f);
extern float fmodf(float l, float r);
extern float sqrtf(float f);
extern float atan2f(float a, float b);
extern float cosf(float f);
extern float sinf(float f);
extern void local_save(void *buf, int len);
extern void local_load(void *buf, int len);
extern unsigned char __heap_base;
#define PAGE_SIZE (1 << 16)
/* -- END WASM -- */

/* -- BEGIN MATH -- */
#define abs(a) (((a) < 0) ? -(a) : (a))

#if 0

  char buf[1 << 8] = {0};
  char *wtr = buf;
  *wtr++ = ' ';
  *wtr++ = '<';
  *wtr++ = '-';
  *wtr++ = ' ';
  fmt_float(wtr, abs(state.hand.x));
                                                
  TextPopup tp = { .msg = buf, .size = 0.03f };
  text_popup_set_in_world(&tp, state.hand.x, state.hand.y);
  text_popup_draw(&tp);
static void fmt_float(char *str, float f) {
  str[0] = ' ';
  str[1] = ' ';
  int i = (int)(f*1000);
  for (int j = 0; j < 4; j++)
    str[2 + 3-1-j] = '0'+(i%10), i/= 10;

  str[0] = str[1];
  str[1] = '.';
}
#endif

static float inline ease_out_quad(float x) { return 1 - (1 - x) * (1 - x); }
static float inline lerp(float v0, float v1, float t) { return (1 - t) * v0 + t * v1; }
static float lerp_rads(float a, float b, float t) {
  float difference = fmodf(b - a, MATH_PI*2.0f),
        distance = fmodf(2.0f * difference, MATH_PI*2.0f) - difference;
  return a + distance * t;
}
static float rads_dist(float a, float b) {
  float difference = fmodf(b - a, MATH_PI*2.0f),
        distance = fmodf(2.0f * difference, MATH_PI*2.0f) - difference;
  return distance;
}

static inline float mag(float x, float y) { return sqrtf(x*x + y*y); }
static inline void norm(float *x, float *y) {
  float m = mag(*x, *y);
  if (m > 0.0f)
    *x /= m,
    *y /= m;
}
#if 0
static inline void reflect(float *_vx, float *_vy, float nx, float ny) {
  float vx = *_vx;
  float vy = *_vy;

  float vdotn = vx*nx + vy*ny;
  *_vx = vx - (2 * vdotn * nx);
  *_vy = vy - (2 * vdotn * ny);
}
#endif

/* https://www.alienryderflex.com/intersect/ */
int line_intersection_free(
  float Ax, float Ay,
  float Bx, float By,
  float Cx, float Cy,
  float Dx, float Dy,
  float *X, float *Y
) {
  float  distAB, theCos, theSin, newX, ABpos ;

  //  Fail if either line is undefined.
  if ((Ax==Bx && Ay==By) || (Cx==Dx && Cy==Dy)) return 0;

  //  (1) Translate the system so that point A is on the origin.
  Bx-=Ax; By-=Ay;
  Cx-=Ax; Cy-=Ay;
  Dx-=Ax; Dy-=Ay;

  //  Discover the length of segment A-B.
  distAB=sqrtf(Bx*Bx+By*By);

  //  (2) Rotate the system so that point B is on the positive X axis.
  theCos=Bx/distAB;
  theSin=By/distAB;
  newX=Cx*theCos+Cy*theSin;
  Cy  =Cy*theCos-Cx*theSin; Cx=newX;
  newX=Dx*theCos+Dy*theSin;
  Dy  =Dy*theCos-Dx*theSin; Dx=newX;

  //  Fail if the lines are parallel.
  if (Cy==Dy) return 0;

  //  (3) Discover the position of the intersection point along line A-B.
  ABpos=Dx+(Cx-Dx)*Dy/(Dy-Cy);

  //  (4) Apply the discovered position to line A-B in the original coordinate system.
  *X=Ax+ABpos*theCos;
  *Y=Ay+ABpos*theSin;

  //  Success.
  return 1;
}

int line_intersection(
  float Ax, float Ay,
  float Bx, float By,
  float Cx, float Cy,
  float Dx, float Dy,
  float *X, float *Y
) {
  float  distAB, theCos, theSin, newX, ABpos ;

  //  Fail if either line segment is zero-length.
  if ((Ax==Bx && Ay==By) || (Cx==Dx && Cy==Dy)) return 0;

  //  Fail if the segments share an end-point.
  if ((Ax==Cx && Ay==Cy) || (Bx==Cx && By==Cy)
  ||  (Ax==Dx && Ay==Dy) || (Bx==Dx && By==Dy)) {
    return 0; }

  //  (1) Translate the system so that point A is on the origin.
  Bx-=Ax; By-=Ay;
  Cx-=Ax; Cy-=Ay;
  Dx-=Ax; Dy-=Ay;

  //  Discover the length of segment A-B.
  distAB=sqrtf(Bx*Bx+By*By);

  //  (2) Rotate the system so that point B is on the positive X axis.
  theCos=Bx/distAB;
  theSin=By/distAB;
  newX=Cx*theCos+Cy*theSin;
  Cy  =Cy*theCos-Cx*theSin; Cx=newX;
  newX=Dx*theCos+Dy*theSin;
  Dy  =Dy*theCos-Dx*theSin; Dx=newX;

  //  Fail if segment C-D doesn't cross line A-B.
  if ((Cy<0. && Dy<0.) || (Cy>=0. && Dy>=0.)) return 0;

  //  (3) Discover the position of the intersection point along line A-B.
  ABpos=Dx+(Cx-Dx)*Dy/(Dy-Cy);

  //  Fail if segment C-D crosses line A-B outside of segment A-B.
  if (ABpos<0. || ABpos>distAB) return 0;

  //  (4) Apply the discovered position to line A-B in the original coordinate system.
  *X=Ax+ABpos*theCos;
  *Y=Ay+ABpos*theSin;

  //  Success.
  return 1;
}

static float point_on_line(float *p_x, float *p_y,
                           float l0_x, float l0_y,
                           float l1_x, float l1_y) {
  float line_len = (l1_x - l0_x)*(l1_x - l0_x) +
                   (l1_y - l0_y)*(l1_y - l0_y);
  float tri_area_x2 = ((*p_x - l0_x) * (l1_x - l0_x)) +
                      ((*p_y - l0_y) * (l1_y - l0_y));

  float U = tri_area_x2/line_len;
  if (U < 0) U = 0;
  if (U > 1) U = 1;
  *p_x = l0_x + (U * (l1_x - l0_x));
  *p_y = l0_y + (U * (l1_y - l0_y));

  return 0.5f;
}
static float point_to_line(float  p_x, float  p_y,
                           float l0_x, float l0_y,
                           float l1_x, float l1_y) {
  float x = p_x;
  float y = p_y;
  point_on_line(&  x, &  y,
                l0_x, l0_y,
                l1_x, l1_y);
  return mag(p_x - x, p_y - y);
}

typedef struct { float x, y; } Pos;
/* tri[0] is joint, rest are bottom */
static void square_to_tri(Pos points[4], Pos tri[3]) {
  /* if we detect the joint,
   * we have 3 points that make up a triangle,
   * so we can find which side isn't filled. */
  int pair_l = 0, pair_r = 1;
  float pair_dist = mag(points[pair_l].x - points[pair_r].x,
                        points[pair_l].y - points[pair_r].y);
  for (int i = 0; i < 4; i++)
    for (int q = 0; q < 4; q++) {
      if (i == q) continue;

      float iq_dist = mag(points[i].x - points[q].x,
                          points[i].y - points[q].y);
      if (iq_dist < pair_dist) {
        pair_dist = iq_dist;
        pair_l = i;
        pair_r = q;
      }
    }

  tri[0].x = lerp(points[pair_l].x, points[pair_r].x, 0.5f);
  tri[0].y = lerp(points[pair_l].y, points[pair_r].y, 0.5f);
  int wtr = 1;
  for (int i = 0; i < 4; i++)
    if (i != pair_l && i != pair_r)
      tri[wtr++] = points[i];
}

/* -- BEGIN RENDER -- */
typedef struct { uint8_t r, g, b, a; } Color;
Color clr_stone = {  80,  80,  80, 255 };
Color clr_board = (Color) { 105,  55,  35, 255 };

static struct {
  int width, height;
  uint8_t *pixels;

  float alpha;

  struct { float x, y; } cam;
} rendr;

static void world_to_screen(float *x, float *y) {
  /* width expressed in terms of height,
   * scale_x tells you how many times longer than it is high */
  float scale_x = rendr.width/rendr.height;
  float scale_y = 1.0f;

  /* first offset all by cam pos */
  *x -= rendr.cam.x;
  *y -= rendr.cam.y;

  /* zoom */
  *x *= scale_x * (rendr.height / 5.0f);
  *y *= scale_y * (rendr.height / 5.0f);

  *x += rendr. width/2;
  *y += rendr.height/2;

  /* positive y is up yall */
  *y = rendr.height - 1 - *y;
}

static void screen_to_world(float *x, float *y) {
  /* width expressed in terms of height,
   * scale_x tells you how many times longer than it is high */
  float scale_x = rendr.width/rendr.height;
  float scale_y = 1.0f;

  /* positive y is up yall */
  *y = rendr.height - 1 - *y;

  *x -= rendr. width/2;
  *y -= rendr.height/2;

  /* zoom */
  *x /= scale_x * (rendr.height / 5.0f);
  *y /= scale_y * (rendr.height / 5.0f);

  /* last offset all by cam pos */
  *x += rendr.cam.x;
  *y += rendr.cam.y;
}

/* -- END RENDER -- */

/* cuz if you\'re board, then you\'re boring */
typedef enum { BoardStage_Uninit, BoardStage_Init, BoardStage_Sticky } BoardStage;
typedef struct Board Board;
struct Board {
  Board *next, *prev;

  BoardStage stage;

  float beg_x, beg_y,
        end_x, end_y;
};
static void find_nearest_board(float x, float y, Board **best_board, float *to_best_board);
static int32_t board_index(Board *);
static Board *board_at_index(int32_t);


typedef enum {
  HandStage_Init,
  HandStage_Grab,
  HandStage_Bring,
  HandStage_Release,
  HandStage_Return,

  HandStage_Frozen,
} HandStage;
typedef struct {
  HandStage stage;

  float x, y;

  uint32_t tick_stage_start,
           tick_stage_end;

  int32_t grabbed_board;
} Hand;
typedef struct {
  float theta, grip_width;
  float x, y;
} HandOut;

struct { float offset, width; } ARC = {
  .offset = (MATH_PI*2)/3,
  .width = MATH_PI*0.9f,
};
float GRIP_WIDTH_GRAB = MATH_PI*0.25f;
float GRIP_WIDTH_RELEASE = MATH_PI*0.8f;

static void hand_out(Hand *hand, uint32_t tick, HandOut *out) {
  out->x = hand->x;
  out->y = hand->y;

  float duration = hand->tick_stage_end - hand->tick_stage_start;
  float elapsed = tick - hand->tick_stage_start;
  float t = (duration > 0) ? (elapsed / duration) : 0;
  if (t < 0) t = 0;
  if (t > 1) t = 1;
  switch (hand->stage) {
    case (HandStage_Init): {
      hand->stage = HandStage_Grab;
      hand->grabbed_board = -1;
      out->theta = ARC.offset + ARC.width/-2;
      hand->tick_stage_start = tick;
      hand->tick_stage_end   = tick + TICK_SECOND;
    } break;

    case (HandStage_Frozen): {
      out->theta = ARC.offset + ARC.width/-2;
      out->grip_width = GRIP_WIDTH_GRAB;
    } break;

    case (HandStage_Grab): {
      // t = ease_out_quad(t);
      out->theta = ARC.offset + ARC.width/-2;
      out->grip_width = lerp_rads(GRIP_WIDTH_RELEASE, GRIP_WIDTH_GRAB, t);

      if (elapsed < (duration + TICK_SECOND*0.2f)) break;
      hand->stage = HandStage_Bring;
      hand->tick_stage_start = tick;
      hand->tick_stage_end   = tick + TICK_SECOND*1.5;

      Board *nearest_board = 0;
      float to_nearest_board = 0;
      find_nearest_board(hand->x + cosf(out->theta),
                         hand->y + sinf(out->theta),
                         &nearest_board, &to_nearest_board);
      if (to_nearest_board < 0.1f)
        hand->grabbed_board = board_index(nearest_board);
    } break;

    case (HandStage_Bring): {
      t = ease_out_quad(t);
      out->theta = ARC.offset + lerp_rads(ARC.width/-2, ARC.width/2, t);
      out->grip_width = GRIP_WIDTH_GRAB;

      if (t < 1) break;
      hand->stage = HandStage_Release;
      hand->tick_stage_start = tick;
      hand->tick_stage_end   = tick + TICK_SECOND*0.5f;
    } break;

    case (HandStage_Release): {
      // t = ease_out_quad(t);
      out->grip_width = lerp_rads(GRIP_WIDTH_GRAB, GRIP_WIDTH_RELEASE, t);
      out->theta = ARC.offset + ARC.width/2;

      if (elapsed < (duration + TICK_SECOND*0.2f)) break;
      if (hand->grabbed_board > -1)
        board_at_index(hand->grabbed_board)->stage = BoardStage_Sticky;
      hand->grabbed_board = -1;
      hand->stage = HandStage_Return;
      hand->tick_stage_start = tick;
      hand->tick_stage_end   = tick + TICK_SECOND;
    } break;

    case (HandStage_Return): {
      t = ease_out_quad(t);
      out->theta = ARC.offset + lerp_rads(ARC.width/2, ARC.width/-2, t);
      out->grip_width = GRIP_WIDTH_RELEASE;

      if (t < 1) break;
      hand->stage = HandStage_Grab;
      hand->tick_stage_start = tick;
      hand->tick_stage_end   = tick + TICK_SECOND*0.5;
    } break;
  }
}

typedef enum {
  TutorialStage_Init,
  TutorialStage_AwaitHandSelect,
  TutorialStage_AwaitHandDragStart,
  TutorialStage_AwaitHandDragEnd,
  TutorialStage_AwaitHouseBuild,
  TutorialStage_AwaitHouseLand,
  TutorialStage_AwaitWoodDispenser,
} TutorialStage;

typedef struct {
  float x, y, size; /* size in pixels = size * rendr.height */
  char *msg;
} TextPopup;

typedef enum {
  WoodDispenserStage_Uninit,
  WoodDispenserStage_Init,
  WoodDispenserStage_Spawning,
  WoodDispenserStage_Holding,
} WoodDispenserStage;
typedef struct {
  WoodDispenserStage stage;
  float x, y;

  uint32_t tick_stage_start, tick_stage_end;
} WoodDispenser;
typedef struct {
  float x, y;
  float beg_x, beg_y,
        end_x, end_y;
} WoodDispenserOut;
static Board *board_alloc(void);
static void text_popup_set_in_world(TextPopup *tp, float x, float y);
static void text_popup_draw(TextPopup *tp);
static void wood_dispenser_out(WoodDispenser *wd, uint32_t tick, WoodDispenserOut *wdo) {
  wdo->x = wd->x;
  wdo->y = wd->y;

  float duration = wd->tick_stage_end - wd->tick_stage_start;
  float elapsed = tick - wd->tick_stage_start;
  float t = (duration > 0) ? (elapsed / duration) : 0;
  if (t < 0) t = 0;
  if (t > 1) t = 1;
  switch (wd->stage) {
    case (WoodDispenserStage_Uninit): break;
    case (WoodDispenserStage_Init): {
      wd->stage = WoodDispenserStage_Spawning;
      wd->tick_stage_start = tick;
      wd->tick_stage_end   = tick + TICK_SECOND;
    } break;
    case (WoodDispenserStage_Spawning): {
      float board_scale = ease_out_quad(t);
      {
        float t = ((double)tick)/TICK_SECOND * 6.0f;
        float dx = cosf(cosf(t)*0.07f);
        float dy = sinf(cosf(t)*0.07f);
        float d = board_scale/2.0f;
        wdo->beg_x = wd->x - dx*d;
        wdo->beg_y = wd->y - dy*d;
        wdo->end_x = wd->x + dx*d;
        wdo->end_y = wd->y + dy*d;
      }

      if (t < 1) break;
      wd->stage = WoodDispenserStage_Holding;
      Board *b = board_alloc();
      if (b) {
        b->stage = BoardStage_Init;
        b->beg_x = wdo->beg_x;
        b->beg_y = wdo->beg_y;
        b->end_x = wdo->end_x;
        b->end_y = wdo->end_y;
        break;
      }
    } break;
    case (WoodDispenserStage_Holding): {
      float dist = 0;
      find_nearest_board(wd->x+0.05f, wd->y-0.05f, 0, &dist);

      
      
      
      
      
      
      

      
      
      

      if (dist > 0.1f)
        wd->stage = WoodDispenserStage_Init;
    } break;
  }

}

typedef struct {
  int init;
  Board boards[2];
} House;

typedef enum {
  Mode_View,
  Mode_HandTweak,
  Mode_HandMove,
  Mode_BuyPreview,
} Mode;
static struct { 
  uint32_t tick, tick_last_click;
  Mode mode;
  /* for Mode_Hand* */ Hand *hand_selected;
  WoodDispenser wood_dispensers[1 << 4];

  TutorialStage tutorial_stage;

  Hand hand;
  Board boards[1 << 5];

  House house;
} state = {0};

static struct {
  float mouse_x, mouse_y;
  double elapsed;
} env;


static void text_popup_set_in_world(TextPopup *tp, float x, float y);
static void text_popup_draw(TextPopup *tp);
static void plot_line(float p0_x, float p0_y,
                      float p1_x, float p1_y, Color c);
static void board_pivot(Board *b, float x, float y, float delta_theta);
static void tutorial(void) {

  float ty = 1.5f;

  float hand_ideal_x = 0.465f;
  float hand_ideal_y = 0.735f+ty;

  switch (state.tutorial_stage) {
    case (TutorialStage_Init): {
#if 0
      local_load(&state.boards, sizeof(state.boards));
      for (int i = 0; i < ARR_LEN(state.boards); i++) {
        if (state.boards[i].next) state.boards[i].next += (long)(state.boards - 1);
        if (state.boards[i].prev) state.boards[i].prev += (long)(state.boards - 1);
      }
      state.house = state.boards;
#else 
      state.hand.x = hand_ideal_x*0.6f;
      state.hand.y = ty;

      state.boards[0] = (Board) {
        .stage = BoardStage_Init,
        .beg_x = 0.8f + 1.2f, .beg_y =        0.6f+ty,
        .end_x =        1.2f, .end_y = 0.8f + 0.6f+ty,
      };

      state.boards[1] = (Board) {
        .stage = BoardStage_Init,
        .beg_x =  0.8f - 1.4f, .beg_y =          1.4f+ty,
        .end_x =       - 1.4f, .end_y = - 0.8f + 1.4f+ty,
      };
#endif

      state.tutorial_stage = TutorialStage_AwaitHandSelect;
    } break;

    case (TutorialStage_AwaitHandSelect): {

      text_popup_draw(&(TextPopup) {
        .msg = "can't reach!",
        .size = 0.12f, .y = 0.1f, .x = 0.2f
      });

      text_popup_draw(&(TextPopup) {
        .msg = "mf got t-rex arm",
        .size = 0.07f, .y = 0.02f, .x = 0.3f
      });

      TextPopup tp = { .msg = " <- double click here to select", .size = 0.03f };
      text_popup_set_in_world(&tp, state.hand.x, state.hand.y);
      text_popup_draw(&tp);

      if (state.mode == Mode_HandTweak) {
        state.tutorial_stage = TutorialStage_AwaitHandDragStart;
      }
    } break;

    case (TutorialStage_AwaitHandDragStart): {
      text_popup_draw(&(TextPopup) {
        .msg = "selected!",
        .size = 0.12f, .y = 0.1f, .x = 0.2f
      });

      text_popup_draw(&(TextPopup) {
        .msg = "now drag 'em so he can reach!",
        .size = 0.07f, .y = 0.02f, .x = 0.1f
      });

      if (state.mode == Mode_HandMove) {
        state.tutorial_stage = TutorialStage_AwaitHandDragEnd;
      }
    } break;

    case (TutorialStage_AwaitHandDragEnd): {
      state.hand_selected->stage = HandStage_Frozen;
      TextPopup tp = { .msg = " <- here maybe?", .size = 0.03f };
      text_popup_set_in_world(&tp, hand_ideal_x, hand_ideal_y);
      text_popup_draw(&tp);

      float dx = state.hand_selected->x - hand_ideal_x;
      float dy = state.hand_selected->y - hand_ideal_y;
      if (mag(dx, dy) > 0.2f) break;

      if (state.mode == Mode_HandMove) break;
      state.hand_selected->stage = HandStage_Init;
      state.hand_selected->x = hand_ideal_x;
      state.hand_selected->y = hand_ideal_y;
      state.tutorial_stage = TutorialStage_AwaitHouseBuild;
    } break;

    case (TutorialStage_AwaitHouseBuild): {
      if (!state.house.init) break;
      state.tutorial_stage = TutorialStage_AwaitHouseLand;
    } break;

    case (TutorialStage_AwaitHouseLand): {
      text_popup_draw(&(TextPopup) {
        .msg = "a home!",
        .size = 0.12f, .y = 0.1f, .x = 0.2f
      });

      text_popup_draw(&(TextPopup) {
        .msg = "we'll be needing more of these ...",
        .size = 0.07f, .y = 0.02f, .x = 0.0f
      });

      Board *bs = state.house.boards;
      // board_pivot(bs+0, cx, cy, 0.1f);
      // board_pivot(bs+1, cx, cy, 0.1f);

      Pos
        tri[3] = {0},
        points[4] = {
          { bs[0].beg_x, bs[0].beg_y },
          { bs[0].end_x, bs[0].end_y },
          { bs[1].beg_x, bs[1].beg_y },
          { bs[1].end_x, bs[1].end_y },
        };
      square_to_tri(points, tri);

      /* center of base of triangle */
      float cx = (tri[1].x + tri[2].x) / 2;
      float cy = (tri[1].y + tri[2].y) / 2;
      float dx = tri[1].x - tri[2].x;
      float dy = tri[1].y - tri[2].y;

      /* intentionally a perpendicular vector (also fuck atan2f's fn signature) */
      float rot_now = atan2f( -dx, dy );
      float ideal_rot = MATH_PI/2.0f;
      float delta_rot = rads_dist(rot_now, ideal_rot);
      board_pivot(bs+0, cx, cy, delta_rot*0.04f);
      board_pivot(bs+1, cx, cy, delta_rot*0.04f);

      bs[0].beg_y -= (cy - (cy*0.95f));
      bs[0].end_y -= (cy - (cy*0.95f));
      bs[1].beg_y -= (cy - (cy*0.95f));
      bs[1].end_y -= (cy - (cy*0.95f));

      if (delta_rot < 0.01f && cy <= 0.01f)
        state.tutorial_stage = TutorialStage_AwaitWoodDispenser;
    } break;

    case (TutorialStage_AwaitWoodDispenser): {
      text_popup_draw(&(TextPopup) {
        .msg = "we'll be needing more of these ...",
        .size = 0.07f, .y = 0.02f, .x = 0.0f
      });

      text_popup_draw(&(TextPopup) {
        .msg = " click here -> ",
        .size = 0.06f, .y = 0.93f, .x = 0.36f
      });
    }
  }
}

static void find_nearest_hand(float x, float y, Hand **best_hand, float *to_best_hand) {
  HandOut hando = {0};

  Hand *hand = &state.hand;
  hand_out(hand, state.tick, &hando);

  float beg_x = state.hand.x;
  float beg_y = state.hand.y;
  float end_x = beg_x + cosf(hando.theta);
  float end_y = beg_y + sinf(hando.theta);
  *to_best_hand = point_to_line(x, y, beg_x, beg_y, end_x, end_y);
  *best_hand = hand;
}

static int32_t board_index(Board *board) { return board - state.boards; }
static Board *board_at_index(int32_t i) {
  return (i >= 0 && i < ARR_LEN(state.boards)) ? state.boards+ i : 0;
}
static void find_nearest_board(float x, float y, Board **best_board, float *to_best_board) {
  if (*to_best_board == 0) *to_best_board = 1e9;

  for (int i = 0; i < ARR_LEN(state.boards); i++) {
    Board *b = state.boards + i;
    if (b->stage == BoardStage_Uninit) continue;

    float to_board = point_to_line(x, y, b->beg_x, b->beg_y, b->end_x, b->end_y);

    if (to_board < *to_best_board) {
      *to_best_board = to_board;
      if (best_board) *best_board = b;
    }
  }
}
static Board *board_alloc(void) {
  for (int i = 0; i < ARR_LEN(state.boards); i++) {
    Board *b = state.boards + i;
    if (b->stage == BoardStage_Uninit) {
      __builtin_memset(b, 0, sizeof(*b));
      return b;
    }
  }
  return 0;
}
static WoodDispenser *wood_dispenser_alloc(void) {
  for (int i = 0; i < ARR_LEN(state.wood_dispensers); i++) {
    WoodDispenser *b = state.wood_dispensers + i;
    if (b->stage == WoodDispenserStage_Uninit) {
      __builtin_memset(b, 0, sizeof(*b));
      return b;
    }
  }
  return 0;
}

WASM_EXPORT void mousemove(float x, float y) {
  screen_to_world(&x, &y);
  env.mouse_x = x;
  env.mouse_y = y;
}

WASM_EXPORT void mouseup(float x, float y) {
  screen_to_world(&x, &y);
  env.mouse_x = x;
  env.mouse_y = y;

  if (state.mode == Mode_HandMove) {
    state.hand_selected->x = x;
    state.hand_selected->y = y;
    state.mode = Mode_HandTweak;
  }
}

float CLICK_DIST_SELECT   = 0.2f;
float CLICK_DIST_UNSELECT = 0.5f;

static void ui(int mousedown);
WASM_EXPORT void mousedown(float x, float y) {
  int double_click = (state.tick - state.tick_last_click) < (TICK_SECOND/3);
  state.tick_last_click = state.tick;

  screen_to_world(&x, &y);
  env.mouse_x = x;
  env.mouse_y = y;

  Hand *nearest_hand = 0;
  float to_nearest_hand = 0;
  find_nearest_hand(x, y, &nearest_hand, &to_nearest_hand);

  switch (state.mode) {
    case (Mode_View): {
      if (to_nearest_hand < CLICK_DIST_SELECT && double_click) {
        state.hand_selected = &state.hand;
        state.mode = Mode_HandTweak;
      }
    } break;

    case (Mode_HandTweak): {
      if (to_nearest_hand > CLICK_DIST_UNSELECT) {
        state.mode = Mode_View;
        state.hand_selected = 0;
      }
      if (to_nearest_hand < CLICK_DIST_SELECT) {
        state.mode = Mode_HandMove;
      }
    } break;

    case (Mode_HandMove): {
      // hi
    } break;

    case (Mode_BuyPreview): {
      // hi (see ui)
    } break;
  }

  ui(1);
}

WASM_EXPORT uint8_t *init(int width, int height) {
  rendr.width = width;
  rendr.height = height;
  rendr.pixels = &__heap_base;
  int mem_needed = (width * height * 4)/PAGE_SIZE;
  int delta = mem_needed - __builtin_wasm_memory_size(0) + 2;

  if (delta > 0) __builtin_wasm_memory_grow(0, delta);

  return rendr.pixels;
}

static void write_pixel(int x, int y, Color src, float a) {
  /* premultiply src */
  // src.r = 255*(1 - (float)(255 - src.r)/255 * a);
  // src.g = 255*(1 - (float)(255 - src.g)/255 * a);
  // src.b = 255*(1 - (float)(255 - src.b)/255 * a);
  a *= rendr.alpha;
  src.r *= a;
  src.g *= a;
  src.b *= a;

  if (x < 0 || x >= rendr. width) return;
  if (y < 0 || y >= rendr.height) return;
  Color *dst = ((Color *)rendr.pixels) + ((y * rendr.width) + x);
  dst->r = src.r + dst->r*(1 - a);
  dst->g = src.g + dst->g*(1 - a);
  dst->b = src.b + dst->b*(1 - a);
  dst->a = 255;
}

#if 0
static void fill_rect(int px, int py) {
  for (int x = px-5; x < px+5; x++)
    for (int y = py-5; y < py+5; y++)
      write_pixel(x, y, (Color) { 255, 255, 9, 255 });
}
#endif

// integer part of x
static int ipart(float f) { return (int)f; }
static int round(float x) { return ipart(x + 0.5); }

// fractional part of x
static float fpart(float x) { return x - ipart(x); }
static float rfpart(float x) { return 1 - fpart(x); }

static void plot_line(float p0_x, float p0_y,
                      float p1_x, float p1_y, Color c) {
  float x0 = p0_x, y0 = p0_y;
  float x1 = p1_x, y1 = p1_y;
  world_to_screen(&x0, &y0);
  world_to_screen(&x1, &y1);
  
  uint8_t steep = abs(y1 - y0) > abs(x1 - x0);

  float tmp;
  #define swap(a, b) tmp = a, a = b, b = tmp
  if   (steep) swap(x0, y0), swap(x1, y1);
  if (x0 > x1) swap(x0, x1), swap(y0, y1);
  #undef swap
  
  float dx = x1 - x0;
  float dy = y1 - y0;

  float gradient = (dx == 0.0f) ? 1.0f : dy / dx;

  // handle first endpoint
  float xend = round(x0);
  float yend = y0 + gradient * (xend - x0);
  float xgap = rfpart(x0 + 0.5);
  float xpxl1 = xend; // this will be used in the main loop
  float ypxl1 = ipart(yend);
  if (steep) {
    write_pixel(ypxl1,   xpxl1, c, rfpart(yend) * xgap);
    write_pixel(ypxl1+1, xpxl1, c,  fpart(yend) * xgap);
  } else {
    write_pixel(xpxl1, ypxl1  , c, rfpart(yend) * xgap);
    write_pixel(xpxl1, ypxl1+1, c,  fpart(yend) * xgap);
  }

  float intery = yend + gradient; // first y-intersection for the main loop

  // handle second endpoint
  xend = round(x1);
  yend = y1 + gradient * (xend - x1);
  xgap = fpart(x1 + 0.5);
  float xpxl2 = xend; // this will be used in the main loop
  float ypxl2 = ipart(yend);
  if (steep) {
    write_pixel(ypxl2  , xpxl2, c, rfpart(yend) * xgap);
    write_pixel(ypxl2+1, xpxl2, c,  fpart(yend) * xgap);
  } else {
    write_pixel(xpxl2, ypxl2,   c, rfpart(yend) * xgap);
    write_pixel(xpxl2, ypxl2+1, c,  fpart(yend) * xgap);
  }

  // main loop
  if (steep) {
    for (int x = xpxl1 + 1; x < xpxl2; x++) {
      write_pixel(ipart(intery)  , x, c, rfpart(intery));
      write_pixel(ipart(intery)+1, x, c,  fpart(intery));
      intery += gradient;
    }
  } else {
    for (int x = xpxl1 + 1; x < xpxl2; x++) {
      write_pixel(x, ipart(intery),   c, rfpart(intery));
      write_pixel(x, ipart(intery)+1, c,  fpart(intery));
      intery += gradient;
    }
  }
}

static void plot_line_thick(
  float p0_x, float p0_y,
  float p1_x, float p1_y,
  Color c, float thickness
) {
  thickness /= 2.0f;

  float perp_x = p1_x - p0_x;
  float perp_y = p1_y - p0_y;
  norm(&perp_x, &perp_y);
  float dx =  perp_y;
  float dy = -perp_x;

  plot_line(p0_x - dx*thickness, p0_y - dy*thickness,
            p1_x - dx*thickness, p1_y - dy*thickness,
            c);
  plot_line(p0_x + dx*thickness, p0_y + dy*thickness,
            p1_x + dx*thickness, p1_y + dy*thickness,
            c);
  plot_line(p1_x + dx*thickness, p1_y + dy*thickness,
            p1_x - dx*thickness, p1_y - dy*thickness,
            c);
  plot_line(p0_x + dx*thickness, p0_y + dy*thickness,
            p0_x - dx*thickness, p0_y - dy*thickness,
            c);
}

static void draw_hand(HandOut *hando, Color body_clr) {
  float FINGER_LEN = 0.2f;
  float theta      = hando->theta;
  float grip_width = hando->grip_width;

  float beg_x = hando->x;
  float beg_y = hando->y;
  float end_x = beg_x + cosf(theta) * (1-FINGER_LEN);
  float end_y = beg_y + sinf(theta) * (1-FINGER_LEN);

  plot_line_thick(beg_x, beg_y, end_x, end_y, body_clr, 0.085f);
  {
    float x = end_x;
    float y = end_y;
    float q;
    q = theta - grip_width/2;
    plot_line_thick(x, y, x + cosf(q)*FINGER_LEN, y + sinf(q)*FINGER_LEN, body_clr, 0.065f);
    q = theta + grip_width/2;
    plot_line_thick(x, y, x + cosf(q)*FINGER_LEN, y + sinf(q)*FINGER_LEN, body_clr, 0.065f);
  }
}

/* x, y in domain 0..8, c is char */
static int fontdata_read(int x, int y, char c) {
  uint8_t bits = fontdata[c*8 + y];
  return (bits >> (7-x)) & 1;
}

static void rcx_char(int px, int py, int scale, char c) {
  for (int y = 0; y < 8*scale; y++) {
    for (int x = 0; x < 8*scale; x++)
      if (fontdata_read(x/scale, y/scale, c))
        write_pixel(px + x, py + y, (Color) { 0, 0, 0, 255 }, 1.0f);
  }
}

static void text_popup_set_in_world(TextPopup *tp, float x, float y) {
  world_to_screen(&x, &y);
  tp->x =     x / rendr. width;
  tp->y = 1 - y / rendr.height;
}
static void text_popup_draw(TextPopup *tp) {
  float ideal_size = rendr.height * tp->size;
  int scale = ideal_size/8;

  float x =      tp->x  * rendr. width;
  float y = (1 - tp->y) * rendr.height;
  y -= 8*(scale*0.55);

  for (char *str = tp->msg; *str; str++) {
    rcx_char(x, y, scale, *str);
    x += 8*scale;
  }
}
/* in normalized (screen) coords */
static float text_popup_width(TextPopup *tp) {
  float ideal_size = rendr.height * tp->size;
  int scale = ideal_size/8;

  /* get in pixels then normalize to screen */
  float ret = 0;
  for (char *str = tp->msg; *str; str++) ret += 8*scale;
  return ret/rendr.width;
}
static int text_popup_hovered(TextPopup *tp) {
  float mx = env.mouse_x;
  float my = env.mouse_y;
  world_to_screen(&mx, &my);
  mx /= rendr.width;
  my /= rendr.height;
  my = 1 - my;

  mx -= tp->x;
  my -= tp->y - tp->size/2;

  if (mx < 0 || mx > text_popup_width(tp)) return 0;
  if (my < 0 || my > tp->size) return 0;
  return 1;
}

static void board_pivot(Board *b, float x, float y, float delta_theta) {
  float old_x, old_y;

  old_x = b->beg_x - x;
  old_y = b->beg_y - y;
  b->beg_x = x + old_x*cosf(delta_theta) - old_y*sinf(delta_theta);
  b->beg_y = y + old_x*sinf(delta_theta) + old_y*cosf(delta_theta);

  old_x = b->end_x - x;
  old_y = b->end_y - y;
  b->end_x = x + old_x*cosf(delta_theta) - old_y*sinf(delta_theta);
  b->end_y = y + old_x*sinf(delta_theta) + old_y*cosf(delta_theta);
}

static void structure_test(Board *head) {
  /* it's probably a house if:
   * - there are two links in the chain
   * - the tips of the boards are close
   * - could do angle, but area seems less brittle */

  if (state.house.init) return;

  if (!head || !head->next || (head->next->next != 0)) return;

  state.house.init = 1;
  state.house.boards[0] = *head;
  state.house.boards[1] = *head->next;

        head->stage = BoardStage_Uninit;
  head->next->stage = BoardStage_Uninit;
}

static void draw_wood_dispenser(WoodDispenserOut *wdo) {
  plot_line_thick(wdo->beg_x, wdo->beg_y,
                  wdo->end_x, wdo->end_y,
                  clr_board, 0.2f); 

  for (int p = 0; p < 3; p++) {
    float flip = ((p%2) ? -1 : 1);

    float px = wdo->x + 0.18f*(p-1);
    float py = wdo->y + flip*0.15f;

    float last_x, last_y;
    for (int i = 0; i < 3+1; i++) {
      float t = (float)i/3.0f;
      float r = t * MATH_TAU + MATH_PI*0.501*-flip;
      float x = px + cosf(r)*0.1,
            y = py + sinf(r)*0.1;
      if (i > 0)
        plot_line(last_x, last_y, x, y, (Color) { 0, 0, 0, 255 });
      last_x = x, last_y = y;
    }
  }
}


WASM_EXPORT void draw(double elapsed) {
  rendr.alpha = 1;
  rendr.cam.y = 1.5f;
  env.elapsed = elapsed;

  // __builtin_memset(rendr.pixels, 235, rendr.width * rendr.height * 4);
  __builtin_memset(rendr.pixels, 255, rendr.width * rendr.height * 4);

  tutorial();

  for (int i = 0; i < ARR_LEN(state.boards); i++) {
    if (state.boards[i].next) state.boards[i].next -= (long)(state.boards - 1);
    if (state.boards[i].prev) state.boards[i].prev -= (long)(state.boards - 1);
  }
  local_save(&state.boards, sizeof(state.boards));
  for (int i = 0; i < ARR_LEN(state.boards); i++) {
    if (state.boards[i].next) state.boards[i].next += (long)(state.boards - 1);
    if (state.boards[i].prev) state.boards[i].prev += (long)(state.boards - 1);
  }

  int last_tick = state.tick++;

  /* boards move with hands */
  {
    Hand *hand = &state.hand;
    HandOut hando = {0};

    hand_out(hand, last_tick, &hando);
    float last_theta = hando.theta;

    hando = (HandOut) {0};

    hand_out(hand, state.tick, &hando);
    float new_theta = hando.theta;

    if (hand->grabbed_board >= 0 && state.boards[hand->grabbed_board].stage) {
      Board *b = state.boards + hand->grabbed_board;

      /* what is this, the name of a fraternity? */
      float delta_theta = new_theta - last_theta;

      Board *head = b; while (head->prev) head = head->prev;
      Board *next = head;
      do {
        board_pivot(next, hand->x, hand->y, delta_theta);
      } while ((next = next->next));
    }
  }

  /* boards that touch, stick */
  for (int i = 0; i < ARR_LEN(state.boards); i++) {
    Board *a = state.boards + i;
    if (a->stage != BoardStage_Sticky) continue;
    a->stage = BoardStage_Init;

    /* quadratic perf goes WEEEE */
    for (int i = 0; i < ARR_LEN(state.boards); i++) {
      Board *b = state.boards + i;
      if (a->stage == BoardStage_Uninit) continue;

      Board *head = b; while (head->prev) head = head->prev;
      Board *next = head;
      do {
        if (next == a) goto SKIP;
      } while ((next = next->next));

      float ox, oy;
      int hit = line_intersection(
        a->beg_x, a->beg_y,
        a->end_x, a->end_y,
        b->beg_x, b->beg_y,
        b->end_x, b->end_y,
        &ox, &oy
      );

      if (hit) {
        if (b->next) b->next->prev = a;

        a->next = b->next;
        b->next = a;

        a->prev = b;

        structure_test(head);
      }
SKIP:
      ;
    }
  }

  /* ... */

  /* ground */
  plot_line(-10.0f, 0.0f,
             10.0f, 0.0f,
            (Color) { 128, 200, 128, 255 });

  {
    HandOut hando = {0};
    hand_out(&state.hand, state.tick, &hando);
    Color body_clr = {  85,  75, 100, 255 };

    Hand *nearest_hand = 0;
    float to_nearest_hand = 0;
    find_nearest_hand(env.mouse_x, env.mouse_y, &nearest_hand, &to_nearest_hand);

    switch (state.mode) {
      case (Mode_View): {
        if (to_nearest_hand < CLICK_DIST_SELECT)
          body_clr = (Color) { 185,  75,  55, 255 };
      } break;

      case (Mode_HandTweak): {
        body_clr = (Color) { 215,  50,   0, 255 };

        {
          float x = hando.x;
          float y = hando.y;
          plot_line(x-0.1, y-0.0, x+0.1, y+0.0, body_clr);
          plot_line(x-0.0, y-0.1, x+0.0, y+0.1, body_clr);

          HandOut arc = hando;
          Color arc_clr = { 205, 140, 205, 255 };
          arc.theta = ARC.offset - ARC.width/2;
          arc.grip_width = GRIP_WIDTH_GRAB;
          draw_hand(&arc, arc_clr);

          arc.theta = ARC.offset + ARC.width/2;
          arc.grip_width = GRIP_WIDTH_RELEASE;
          draw_hand(&arc, arc_clr);
        }

        if (to_nearest_hand > CLICK_DIST_UNSELECT)
          body_clr = (Color) { 180,  85,   0, 255 };
        if (to_nearest_hand < CLICK_DIST_SELECT)
          body_clr = (Color) { 255,   0,   0, 255 };
      } break;

      case (Mode_HandMove): {
        body_clr = (Color) {  80, 195,  80, 255 };

        HandOut arc = hando;
        arc.x = env.mouse_x;
        arc.y = env.mouse_y;

        Color arc_clr = { 140, 205, 185, 255 };
        arc.theta = ARC.offset - ARC.width/2;
        arc.grip_width = GRIP_WIDTH_GRAB;
        draw_hand(&arc, arc_clr);

        arc.theta = ARC.offset + ARC.width/2;
        arc.grip_width = GRIP_WIDTH_RELEASE;
        draw_hand(&arc, arc_clr);
      } break;

      case (Mode_BuyPreview): {
        WoodDispenserOut wdo = {0};
        wood_dispenser_out(
          &(WoodDispenser) { .x = env.mouse_x,
                             .y = env.mouse_y,
                             .tick_stage_start = state.tick-TICK_SECOND,
                             .tick_stage_end   = state.tick+1,
                             .stage = WoodDispenserStage_Spawning },
          state.tick,
          &wdo
        );
        rendr.alpha = 0.3f;
        draw_wood_dispenser(&wdo);
        rendr.alpha = 1.0f;
      } break;
    }

    draw_hand(&hando, body_clr);
  }

  if (state.mode == Mode_HandMove) {
    HandOut ghost = {0};
    hand_out(state.hand_selected, state.tick, &ghost);
    ghost.x = env.mouse_x;
    ghost.y = env.mouse_y;
    draw_hand(&ghost, (Color) { 120, 185, 120, 255 });
  }

  for (int i = 0; i < ARR_LEN(state.wood_dispensers); i++) {
    WoodDispenser *wd = state.wood_dispensers + i;
    if (wd->stage == WoodDispenserStage_Uninit) continue;

    WoodDispenserOut wdo = {0};
    wood_dispenser_out(wd, state.tick, &wdo);
    draw_wood_dispenser(&wdo);
  }

  for (int i = 0; i < ARR_LEN(state.boards); i++) {
    Board *b = state.boards + i;
    if (b->stage == BoardStage_Uninit) continue;
    plot_line_thick(b->beg_x, b->beg_y,
                    b->end_x, b->end_y,
                    clr_board, 0.18f);
  }

  {
    float x = env.mouse_x;
    float y = env.mouse_y;
    Color red = { 255, 0, 0, 255 };
    plot_line(x-0.1, y-0.0, x+0.1, y+0.0, red);
    plot_line(x-0.0, y-0.1, x+0.0, y+0.1, red);
  }

  /* decorate house */
  if (state.house.init) {
    Board *b = state.house.boards + 0;
    Board *o = state.house.boards + 1;

    plot_line_thick(b->beg_x, b->beg_y,
                    b->end_x, b->end_y,
                    clr_board, 0.18f);
    plot_line_thick(o->beg_x, o->beg_y,
                    o->end_x, o->end_y,
                    clr_board, 0.18f);

    Pos
      tri[3] = {0},
      points[4] = {
        { b->beg_x, b->beg_y },
        { b->end_x, b->end_y },
        { o->beg_x, o->beg_y },
        { o->end_x, o->end_y },
      };

    square_to_tri(points, tri);
    // Color red = { 255, 0, 0, 255 };
    // plot_line(tri[0].x, tri[0].y,
    //           tri[1].x, tri[1].y, red);
    // 
    // plot_line(tri[1].x, tri[1].y,
    //           tri[2].x, tri[2].y, red);

    /* from one unconnected end to the other */
    float dx = tri[1].x - tri[2].x;
    float dy = tri[1].y - tri[2].y;

    /* from the middle of the unconnected end to the joint */
    float px = lerp(tri[1].x, tri[2].x, 0.5f);
    float py = lerp(tri[1].y, tri[2].y, 0.5f);
    float pnx = tri[0].x - px;
    float pny = tri[0].y - py;
    float pmag = mag(pnx, pny);
    norm(&pnx, &pny);

    /* how long is the bottom of the house? */
    float dmag = mag(dx, dy);
    norm(&dx, &dy);

    /* find the perpendicular vector that points toward the joint */
    float nx, ny;
    nx = dy, ny = -dx;
    if ((nx*pnx+ny*pny) < 0.0f)
      nx = -dy, ny = dx;

    {
      /* chimney */
      float x = lerp(tri[2].x, tri[0].x, 0.4f);
      float y = lerp(tri[2].y, tri[0].y, 0.4f);

      plot_line_thick(x + nx*0.05f, y + ny*0.05f,
                      x + nx*0.40f, y + ny*0.40f, clr_stone, 0.16f);
      plot_line_thick(x + nx*0.35f, y + ny*0.35f,
                      x + nx*0.45f, y + ny*0.45f, clr_stone, 0.23f);

      /* smoke */
      int pi_MAX = 5;
      for (int pi_i = 0; pi_i < pi_MAX; pi_i++) {
        float pi = (float)pi_i + fmodf(elapsed, 1.0f);

        float piq = ease_out_quad(pi/(float)pi_MAX)*pi;
        float px = x + nx*(0.65f + piq*0.10f) + dx*pi*0.08f;
        float py = y + ny*(0.65f + piq*0.10f) + dy*pi*0.08f;
        float radius = 0.1f + pi*0.025f;

        uint8_t w = 170 + 12*pi;
        if (pi < 1) w += (255 - w) * (1 - pi);
        if (pi > (pi_MAX-1)) w += (255 - w) * (pi - (pi_MAX-1));
        Color clr_smoke = { w, w, w, 255 };

        float last_x, last_y;
        for (int i = 0; i < 3+1; i++) {
          float t = (float)i/3.0f;
          float r = t * MATH_TAU + elapsed*0.05f + pi;
          float x = px + cosf(r)*radius,
                y = py + sinf(r)*radius;
          if (i > 0)
            plot_line(last_x, last_y, x, y, clr_smoke);
          last_x = x, last_y = y;
        }
      }
    }

    /* draw horizontal lines from 0<->1 side to 0<->2 side */
    /* (remember, 0 is joint) */
    for (float t = 0; t < (pmag - 0.2f); t += 0.165f) {
      float x1 = lerp(tri[1].x, tri[0].x, t);
      float y1 = lerp(tri[1].y, tri[0].y, t);

      float x2 = lerp(tri[2].x, tri[0].x, t);
      float y2 = lerp(tri[2].y, tri[0].y, t);

      float pad = 0.230f;
      float dist = mag(x1 - x2, y1 - y2);
      if (dist < pad*2.0f) continue;

      float d12x = x2 - x1;
      float d12y = y2 - y1;
      norm(&d12x, &d12y);
      x1 += d12x*pad/2;
      y1 += d12y*pad/2;
      x2 -= d12x*pad;
      y2 -= d12y*pad;
      
      plot_line(x1, y1, x2, y2, clr_stone);
    }

    /* draw vertical lines up from the floor until you hit the roof */
    for (float t = 0; t < dmag; t += 0.205f) {
      float x = tri[2].x + dx*t;
      float y = tri[2].y + dy*t;
      float ox, oy;

      float bx = x - nx*0.3f;
      float by = y - ny*0.3f;

      float d = 50.0f;
      line_intersection_free(
        tri[0].x, tri[0].y,
        tri[1].x, tri[1].y,
              bx,       by,
        nx*d + x, ny*d + y,
        &ox, &oy
      );
      /* clamp d by hit */
      float imag = mag(x - ox, y - oy);
      if (imag < d) d = imag;

      line_intersection_free(
        tri[0].x, tri[0].y,
        tri[2].x, tri[2].y,
              bx,       by,
        nx*d + x, ny*d + y,
        &ox, &oy
      );
      /* clamp d by hit */
      imag = mag(x - ox, y - oy);
      if (imag < d) d = imag;

      d -= 0.325f;
      if (d < 0) d = 0;
      plot_line(x, y, x + nx*d, y + ny*d, clr_stone);
    }
  }

  ui(0);

}

static void ui(int mousedown) {
  TextPopup tp = {
    .msg = "wood dispenser",
    .size = 0.05f,
    .y = 0.95f,
    .x = 0.78f
  };
  text_popup_draw(&tp);
  if (text_popup_hovered(&tp)) {
    float x0 =rendr. width * (tp.x);
    float y0 =rendr.height * (1 - tp.y + tp.size/3);
    screen_to_world(&x0, &y0);
    float x1 =rendr. width * (tp.x + text_popup_width(&tp));
    float y1 =rendr.height * (1 - tp.y + tp.size/3);
    screen_to_world(&x1, &y1);

    plot_line(x0, y0, x1, y1, (Color) { 0, 0, 0, 255 });

    if (mousedown) state.mode = Mode_BuyPreview;
  }

  if (!text_popup_hovered(&tp) && mousedown && state.mode == Mode_BuyPreview) {

    WoodDispenser *wd = wood_dispenser_alloc();
    if (!wd) return;

    wd->stage = WoodDispenserStage_Init;
    wd->x = env.mouse_x;
    wd->y = env.mouse_y;

    state.mode = Mode_View;
  }
}

