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

/* -- BEGIN SERIALIZATION -- */
typedef enum {
  LoVe_Save,
  LoVe_Load,
} LoVe;
typedef struct {
  LoVe lo_ve;
  uint8_t *bytes;
  uint32_t bytes_len;
  uint8_t *cursor;

  int oom;
} LoVeEnv;
static void save_byte(LoVeEnv *lenv, uint8_t u8) {
  if ((lenv->cursor - lenv->bytes) >= lenv->bytes_len)
    lenv->oom = 1;
  else
    *lenv->cursor++ = u8;
}
static uint8_t load_byte(LoVeEnv *lenv) {
  if ((lenv->cursor - lenv->bytes) >= lenv->bytes_len) {
    lenv->oom = 1;
    return 0;
  } else
    return *lenv->cursor++;
}
static void lo_ve_uint32_t(LoVeEnv *lenv, uint32_t *u32) {
  if (lenv->lo_ve == LoVe_Save) {
    save_byte(lenv, (*u32 >>  0) & 0xff);
    save_byte(lenv, (*u32 >>  8) & 0xff);
    save_byte(lenv, (*u32 >> 16) & 0xff);
    save_byte(lenv, (*u32 >> 24) & 0xff);
  }
  if (lenv->lo_ve == LoVe_Load) {
    *u32 = ((uint32_t)load_byte(lenv) <<  0) |
           ((uint32_t)load_byte(lenv) <<  8) |
           ((uint32_t)load_byte(lenv) << 16) |
           ((uint32_t)load_byte(lenv) << 24);
  }
}
static void lo_ve_float(LoVeEnv *lenv, float *f32) {
  uint32_t u32 = 0;
  __builtin_memcpy(&u32, f32, sizeof(u32));
  lo_ve_uint32_t(lenv, &u32);
  __builtin_memcpy(f32, &u32, sizeof(u32));
  // print(*f32);
}
/* -- END SERIALIZATION -- */

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

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in the floats i_x and i_y.
int line_intersection(
  float p0_x, float p0_y,
  float p1_x, float p1_y, 
  float p2_x, float p2_y,
  float p3_x, float p3_y,
  float *i_x, float *i_y
) {
  float s1_x, s1_y, s2_x, s2_y;
  s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
  s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

  float s, t;
  s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
  t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

  if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
    // Collision detected
    if (i_x != 0)
        *i_x = p0_x + (t * s1_x);
    if (i_y != 0)
        *i_y = p0_y + (t * s1_y);
    return 1;
  }

  return 0; // No collision
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
static void quad_to_tri(Pos points[4], Pos tri[3]) {
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
Color clr_hand = { 85,  75, 100, 255 };
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
typedef enum {
  BoardStage_Uninit,
  BoardStage_Init,
  BoardStage_Sticky,
  BoardStage_Fading,
} BoardStage;
typedef struct Board Board;
struct Board {
  Board *next, *head;

  BoardStage stage;

  uint32_t tick_stage_start,
           tick_stage_end;

  float beg_x, beg_y,
        end_x, end_y;
};
static void find_nearest_board(float x, float y, Board **best_board, float *to_best_board);
static int32_t board_index(Board *);
static Board *board_at_index(int32_t);


typedef enum {
  HandStage_Uninit,
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

  float angle_grab, angle_release;

  uint32_t tick_stage_start,
           tick_stage_end;

  int32_t grabbed_board;
} Hand;
static void lo_ve_hand(LoVeEnv *lenv, Hand *hand) {
  lo_ve_uint32_t(lenv, &hand->stage);
  lo_ve_float(lenv, &hand->x);
  lo_ve_float(lenv, &hand->y);
  lo_ve_float(lenv, &hand->angle_grab);
  lo_ve_float(lenv, &hand->angle_release);
  lo_ve_uint32_t(lenv, &hand->tick_stage_start);
  lo_ve_uint32_t(lenv, &hand->tick_stage_end);
  lo_ve_uint32_t(lenv, (uint32_t *)&hand->grabbed_board);
}

typedef struct {
  float theta, grip_width;
  float x, y;
  float angle_grab, angle_release;
} HandOut;

float GRIP_WIDTH_GRAB = MATH_PI*0.25f;
float GRIP_WIDTH_RELEASE = MATH_PI*0.8f;

static void hand_out(Hand *hand, uint32_t tick, HandOut *out) {
  out->x = hand->x;
  out->y = hand->y;
  out->angle_grab = hand->angle_grab;
  out->angle_release = hand->angle_release;

  float duration = hand->tick_stage_end - hand->tick_stage_start;
  float elapsed = tick - hand->tick_stage_start;
  float t = (duration > 0) ? (elapsed / duration) : 0;
  if (t < 0) t = 0;
  if (t > 1) t = 1;
  switch (hand->stage) {
    case (HandStage_Uninit): {} break;

    case (HandStage_Init): {
      hand->stage = HandStage_Grab;
      hand->grabbed_board = -1;
      out->theta = hand->angle_grab;
      hand->tick_stage_start = tick;
      hand->tick_stage_end   = tick + TICK_SECOND;

      float arc_offset = (MATH_PI*2)/3;
      float arc_width = MATH_PI*0.9f;
      hand->angle_grab    = arc_offset + arc_width/-2;
      hand->angle_release = arc_offset + arc_width/ 2;
    } break;

    case (HandStage_Frozen): {
      out->theta = hand->angle_grab;
      out->grip_width = GRIP_WIDTH_GRAB;
    } break;

    case (HandStage_Grab): {
      // t = ease_out_quad(t);
      out->theta = hand->angle_grab;
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
      out->theta = lerp_rads(hand->angle_grab, hand->angle_release, t);
      out->grip_width = GRIP_WIDTH_GRAB;

      if (t < 1) break;
      hand->stage = HandStage_Release;
      hand->tick_stage_start = tick;
      hand->tick_stage_end   = tick + TICK_SECOND*0.5f;
    } break;

    case (HandStage_Release): {
      // t = ease_out_quad(t);
      out->grip_width = lerp_rads(GRIP_WIDTH_GRAB, GRIP_WIDTH_RELEASE, t);
      out->theta = hand->angle_release;

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
      out->theta = lerp_rads(hand->angle_release, hand->angle_grab, t);
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
static void lo_ve_wood_dispenser(LoVeEnv *lenv, WoodDispenser *wd) {
  lo_ve_uint32_t(lenv, &wd->stage);
  lo_ve_float(lenv, &wd->x);
  lo_ve_float(lenv, &wd->y);
  lo_ve_uint32_t(lenv, &wd->tick_stage_start);
  lo_ve_uint32_t(lenv, &wd->tick_stage_end);
}

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

typedef enum {
  HouseStage_Uninit,
  HouseStage_Init,
  HouseStage_Falling,
  HouseStage_Settled,
} HouseStage;
typedef struct {
  HouseStage stage;
  Board boards[2];
} House;

typedef enum {
  ToolKind_WoodDispenser,
  ToolKind_Hand,
} ToolKind;
/*
typedef struct {
  ToolKind kind;
  union {
    Hand *hand;
    WoodDispenser *wood_dispenser;
  } ptr;
} Tool;
*/

typedef enum {
  Mode_View,
  Mode_HandView,
  Mode_HandMoveAngleGrab,
  Mode_HandMoveAngleRelease,
  Mode_HandMove,
  Mode_WoodDispenserView,
  Mode_WoodDispenserMove,
  Mode_BuyPreview,
} Mode;
typedef struct {
  uint32_t tick, tick_last_click;
  TutorialStage tutorial_stage;

  Mode mode;
  /* for Mode_Hand*           */ Hand *hand_selected;
  /* for Mode_WoodDispenser*  */ WoodDispenser *wood_dispenser_selected;
  /* for Mode_BuyPreview      */ ToolKind preview_tool_kind;

  Hand hands[1 << 4];
  WoodDispenser wood_dispensers[1 << 4];
  Board boards[1 << 5];
  House houses[1 << 4];

} State;
static void lo_ve_state(LoVeEnv *lenv, State *state) {
  lo_ve_uint32_t(lenv, &state->tick);
  lo_ve_uint32_t(lenv, &state->tick_last_click);
  lo_ve_uint32_t(lenv, &state->tutorial_stage);
  lo_ve_uint32_t(lenv, &state->mode);

  uint32_t hand_selected_i = state->hand_selected - state->hands;
  lo_ve_uint32_t(lenv, &hand_selected_i);
  state->hand_selected = state->hands + hand_selected_i;

  uint32_t wood_dispenser_selected_i =
    state->wood_dispenser_selected - state->wood_dispensers;
  lo_ve_uint32_t(lenv, &wood_dispenser_selected_i);
  state->wood_dispenser_selected = state->wood_dispensers + wood_dispenser_selected_i;

  lo_ve_uint32_t(lenv, &state->preview_tool_kind);

  uint32_t n_hands = ARR_LEN(state->hands);
  lo_ve_uint32_t(lenv, &n_hands);
  for (int i = 0; i < n_hands; i++)
    lo_ve_hand(lenv, state->hands + i);

  uint32_t n_wood_dispenser = ARR_LEN(state->wood_dispensers);
  lo_ve_uint32_t(lenv, &n_wood_dispenser);
  for (int i = 0; i < n_wood_dispenser; i++)
    lo_ve_wood_dispenser(lenv, state->wood_dispensers + i);
}

static State state = {0};
static void state_save(void) {
  uint8_t buf[sizeof(state)] = {0};
  LoVeEnv lenv = {
    .lo_ve = LoVe_Save,
    .bytes     =        buf,
    .bytes_len = sizeof(buf),
    .cursor    =        buf,
  };
  lo_ve_state(&lenv, &state);
  local_save(buf, sizeof(buf));
}
static void state_load(void) {
  uint8_t buf[sizeof(state)] = {0};
  LoVeEnv lenv = {
    .lo_ve = LoVe_Load,
    .bytes     =        buf,
    .bytes_len = sizeof(buf),
    .cursor    =        buf,
  };
  local_load(buf, sizeof(buf));
  lo_ve_state(&lenv, &state);
}

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
      state.hands[0].stage = HandStage_Init;
      state.hands[0].x = hand_ideal_x*0.6f;
      state.hands[0].y = ty;

      Board *l = board_alloc();
      l->stage = BoardStage_Init;
      l->beg_x = 0.8f + 1.2f; l->beg_y =        0.6f+ty;
      l->end_x =        1.2f; l->end_y = 0.8f + 0.6f+ty;

      Board *r = board_alloc();
      r->stage = BoardStage_Init;
      r->beg_x =  0.8f - 1.4f; r->beg_y =          1.4f+ty;
      r->end_x =       - 1.4f; r->end_y = - 0.8f + 1.4f+ty;
#endif

      state.tutorial_stage = TutorialStage_AwaitHandSelect;
    } break;

    case (TutorialStage_AwaitHandSelect): {

      text_popup_draw(&(TextPopup) {
        .msg = "can't reach!",
        .size = 0.095f, .y = 0.13f, .x = 0.2f
      });

      text_popup_draw(&(TextPopup) {
        .msg = "mf got t-rex arm",
        .size = 0.063f, .y = 0.04f, .x = 0.24f
      });

      TextPopup tp = { .msg = " <- double click here to select", .size = 0.025f };
      text_popup_set_in_world(&tp, state.hands[0].x, state.hands[0].y);
      text_popup_draw(&tp);

      if (state.mode == Mode_HandView) {
        state.tutorial_stage = TutorialStage_AwaitHandDragStart;
      }
    } break;

    case (TutorialStage_AwaitHandDragStart): {
      text_popup_draw(&(TextPopup) {
        .msg = "selected!",
        .size = 0.095f, .y = 0.13f, .x = 0.3f
      });

      text_popup_draw(&(TextPopup) {
        .msg = "now drag 'em so he can reach!",
        .size = 0.053f, .y = 0.04f, .x = 0.10f
      });

      if (state.mode == Mode_HandMove) {
        state.tutorial_stage = TutorialStage_AwaitHandDragEnd;
      }
    } break;

    case (TutorialStage_AwaitHandDragEnd): {
      state.hands[0].stage = HandStage_Frozen;
      TextPopup tp = { .msg = " <- here maybe?", .size = 0.03f };
      text_popup_set_in_world(&tp, hand_ideal_x, hand_ideal_y);
      text_popup_draw(&tp);

      float dx = state.hands[0].x - hand_ideal_x;
      float dy = state.hands[0].y - hand_ideal_y;
      if (mag(dx, dy) > 0.2f) break;

      if (state.mode == Mode_HandMove) break;
      state.hands[0].stage = HandStage_Init;
      state.hands[0].x = hand_ideal_x;
      state.hands[0].y = hand_ideal_y;
      state.tutorial_stage = TutorialStage_AwaitHouseBuild;
    } break;

    case (TutorialStage_AwaitHouseBuild): {
      if (state.houses[0].stage == HouseStage_Uninit) break;
      state.tutorial_stage = TutorialStage_AwaitHouseLand;
    } break;

    case (TutorialStage_AwaitHouseLand): {
      text_popup_draw(&(TextPopup) {
        .msg = "a home!",
        .size = 0.095f, .y = 0.13f, .x = 0.3f
      });

      text_popup_draw(&(TextPopup) {
        .msg = "we'll be needing more of these ...",
        .size = 0.041f, .y = 0.04f, .x = 0.14f
      });

      if (state.houses[0].stage == HouseStage_Settled) 
        state.tutorial_stage = TutorialStage_AwaitWoodDispenser;
    } break;

    case (TutorialStage_AwaitWoodDispenser): {
      text_popup_draw(&(TextPopup) {
        .msg = "we'll be needing more of these ...",
        .size = 0.041f, .y = 0.04f, .x = 0.14f
      });

      text_popup_draw(&(TextPopup) {
        .msg = " click here -> ",
        .size = 0.06f, .y = 0.952f, .x = 0.31f
      });
    }
  }
}

static House *house_alloc(void) {
  for (int i = 0; i < ARR_LEN(state.houses); i++) {
    House *b = state.houses + i;
    if (b->stage == HouseStage_Uninit) {
      __builtin_memset(b, 0, sizeof(*b));
      return b;
    }
  }
  return 0;
}
static Hand *hand_alloc(void) {
  for (int i = 0; i < ARR_LEN(state.hands); i++) {
    Hand *b = state.hands + i;
    if (b->stage == HandStage_Uninit) {
      __builtin_memset(b, 0, sizeof(*b));
      return b;
    }
  }
  return 0;
}
static void find_nearest_hand(float x, float y, Hand **best_hand, float *to_best_hand) {
  if (*to_best_hand == 0) *to_best_hand = 1e9;

  for (int i = 0; i < ARR_LEN(state.hands); i++) {
    Hand *hand = state.hands + i;
    if (hand->stage == HandStage_Uninit) continue;

    HandOut hando = {0};
    hand_out(hand, state.tick, &hando);

    float beg_x = hando.x;
    float beg_y = hando.y;
    float end_x = beg_x + cosf(hando.theta);
    float end_y = beg_y + sinf(hando.theta);
    float to_hand = point_to_line(x, y, beg_x, beg_y, end_x, end_y);

    if (to_hand < *to_best_hand) {
      *to_best_hand = to_hand;
      if (best_hand) *best_hand = hand;
    }
  }
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
      b->head = b;
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
static void find_nearest_wood_dispenser(float x, float y, WoodDispenser **best_wood_dispenser, float *to_best_wood_dispenser) {
  if (*to_best_wood_dispenser == 0) *to_best_wood_dispenser = 1e9;

  for (int i = 0; i < ARR_LEN(state.wood_dispensers); i++) {
    WoodDispenser *b = state.wood_dispensers + i;
    if (b->stage == WoodDispenserStage_Uninit) continue;

    float to_wood_dispenser = mag(x - b->x, y - b->y);;

    if (to_wood_dispenser < *to_best_wood_dispenser) {
      *to_best_wood_dispenser = to_wood_dispenser;
      if (best_wood_dispenser) *best_wood_dispenser = b;
    }
  }
}

WASM_EXPORT void mousemove(float x, float y) {
  state_load();

  screen_to_world(&x, &y);
  env.mouse_x = x;
  env.mouse_y = y;

  state_save();
}

WASM_EXPORT void mouseup(float x, float y) {
  state_load();
  
  screen_to_world(&x, &y);
  env.mouse_x = x;
  env.mouse_y = y;

  WoodDispenser *wd = state.wood_dispenser_selected;
  if (state.mode == Mode_WoodDispenserMove) {
    wd->x = x;
    wd->y = y;
    state.mode = Mode_WoodDispenserView;
  }

  Hand *hand = state.hand_selected;
  if (state.mode == Mode_HandMove) {
    hand->x = x;
    hand->y = y;
    state.mode = Mode_HandView;
  }
  if (state.mode == Mode_HandMoveAngleGrab) {
    hand->angle_grab = atan2f(y - hand->y, x - hand->x);
    hand->stage = HandStage_Grab;
    hand->tick_stage_start = state.tick;
    hand->tick_stage_end   = state.tick + TICK_SECOND;
    state.mode = Mode_HandView;
  }
  if (state.mode == Mode_HandMoveAngleRelease) {
    hand->angle_release = atan2f(y - hand->y, x - hand->x);
    hand->stage = HandStage_Release;
    hand->tick_stage_start = state.tick;
    hand->tick_stage_end   = state.tick + TICK_SECOND;
    state.mode = Mode_HandView;
  }

  state_save();
}

float CLICK_DIST_SELECT   = 0.2f;
float CLICK_DIST_UNSELECT = 0.5f;
 
typedef enum {
  UiEvent_Mousedown,
  UiEvent_DoubleClick,
  UiEvent_Render,
} UiEvent;
static void ui(UiEvent event);
WASM_EXPORT void mousedown(float x, float y) {
  state_load();

  int double_click = (state.tick - state.tick_last_click) < (TICK_SECOND/3);
  state.tick_last_click = state.tick;

  screen_to_world(&x, &y);
  env.mouse_x = x;
  env.mouse_y = y;

  ui(double_click ? UiEvent_DoubleClick : UiEvent_Mousedown);

  state_save();
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

#if 0
static void rcx_char(int px, int py, int scale, char c) {
  for (int y = 0; y < 8*scale; y++) {
    for (int x = 0; x < 8*scale; x++)
      if (fontdata_read(x/scale, y/scale, c))
        write_pixel(px + x, py + y, (Color) { 0, 0, 0, 255 }, 1.0f);
  }
}
#endif

static void text_popup_set_in_world(TextPopup *tp, float x, float y) {
  world_to_screen(&x, &y);
  tp->x =     x / rendr. width;
  tp->y = 1 - y / rendr.height;
}
static void text_popup_draw(TextPopup *tp) {
#if 0
  /* scale by integer factor */
  int ideal_size = rendr.height * tp->size;
  int scale = ideal_size/13.5;
  if (scale < 1) scale = 1;

  float x =      tp->x  * rendr. width;
  float y = (1 - tp->y) * rendr.height;
  y -= 8*(scale*0.55);

  for (char *str = tp->msg; *str; str++) {
    rcx_char(x, y, scale, *str);
    x += 8*scale;
  }
#endif
  int height = rendr.height * tp->size;

  int len = 0;
  for (char *str = tp->msg; *str; str++) len++;
  int width = len*height; /* monospace */

  int nx =      tp->x  * rendr. width;
  int ny = (1 - tp->y) * rendr.height;
  ny -= height*0.55f;

  for (int x = 0; x < width; x++)
    for (int y = 0; y < height; y++) {
      // write_pixel(nx + x, ny + y, (Color) { 0, 0, 0, 255 }, 1.0f);
      // continue;

      int char_i = x / height; /* monospace */

      /* premultiplying by eight instead of using the FPU because im a fuckin chad */
      int char_x = ((x % height) * 8) / height;
      int char_y = ((y * 8) / height);

      if (fontdata_read(char_x, char_y, tp->msg[char_i]))
        write_pixel(nx + x, ny + y, (Color) { 0, 0, 0, 255 }, 1.0f);
    }
}
/* in normalized (screen) coords */
static float text_popup_width(TextPopup *tp) {
#if 0
  old integer scale code 
  float ideal_size = rendr.height * tp->size;
  int scale = ideal_size/8;

  /* get in pixels then normalize to screen */
  float ret = 0;
  for (char *str = tp->msg; *str; str++) ret += 8*scale;
  return ret/rendr.width;
#endif 
  int len = 0;
  for (char *str = tp->msg; *str; str++) len++;
  return len*tp->size * 0.5f; /* monospace */
}
static void text_popup_underline(TextPopup *tp) {
  float x0 =rendr. width * (tp->x);
  float y0 =rendr.height * (1 - tp->y + tp->size/3);
  screen_to_world(&x0, &y0);
  float x1 =rendr. width * (tp->x + text_popup_width(tp));
  float y1 =rendr.height * (1 - tp->y + tp->size/3);
  screen_to_world(&x1, &y1);

  plot_line(x0, y0, x1, y1, (Color) { 0, 0, 0, 255 });
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
   * - the "joint" is tight? (tip touches tip, not center of board)
   * - could do angle, but area seems less brittle? */

  if (!head || !head->next || (head->next->next != 0)) return;

  /* area check */
  {
    Board *b = head;
    Board *o = head->next;
    Pos
      tri[3] = {0},
      points[4] = {
        { b->beg_x, b->beg_y },
        { b->end_x, b->end_y },
        { o->beg_x, o->beg_y },
        { o->end_x, o->end_y },
      };
    quad_to_tri(points, tri);

    float cx = (tri[1].x + tri[2].x) / 2;
    float cy = (tri[1].y + tri[2].y) / 2;
    float area = mag(tri[1].x - tri[2].x, tri[1].y - tri[2].y) *
                 mag(tri[0].x -       cx, tri[0].y -       cy) * 0.5f;
    print(area);

    if (area < 0.35f) return;
  }

  House *house = house_alloc();
  if (house) {
    house->stage = HouseStage_Init;
    house->boards[0] = *head;
    house->boards[1] = *head->next;

          head->stage = BoardStage_Uninit;
    head->next->stage = BoardStage_Uninit;
  }
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
  state_load();

  rendr.alpha = 1;
  rendr.cam.y = 1.5f;
  env.elapsed = elapsed;

  // __builtin_memset(rendr.pixels, 235, rendr.width * rendr.height * 4);
  __builtin_memset(rendr.pixels, 255, rendr.width * rendr.height * 4);

  // tutorial();

  int last_tick = state.tick++;

  /* boards move with hands */
  for (int i = 0; i < ARR_LEN(state.hands); i++) {
    Hand *hand = state.hands + i;
    if (hand->stage == HandStage_Uninit) continue;

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

      Board *next = b->head;
      do {
        /* don't die in my arms ... die after that lol */
        /* no dying while i'm dragging you around!!! */
        next->stage = BoardStage_Fading;
        next->tick_stage_start = state.tick;
        next->tick_stage_end   = state.tick + TICK_SECOND;

        board_pivot(next, hand->x, hand->y, delta_theta);
      } while ((next = next->next));
    }
  }

  for (int i = 0; i < ARR_LEN(state.boards); i++) {
    Board *a = state.boards + i;

    if (a->stage == BoardStage_Fading && state.tick > a->tick_stage_end)
      a->stage = BoardStage_Uninit;

    /* boards that touch, stick */
    if (a->stage != BoardStage_Sticky) continue;
    a->stage = BoardStage_Fading;
    a->tick_stage_start = state.tick;
    a->tick_stage_end   = state.tick + TICK_SECOND;

    /* quadratic perf goes WEEEE */
    for (int i = 0; i < ARR_LEN(state.boards); i++) {
      Board *b = state.boards + i;
      if (a->stage == BoardStage_Uninit) continue;

      Board *next = b->head;
      do {
        if (next == a) goto SKIP;
      } while (next->next && (next = next->next));

      float ox, oy;
      int hit = line_intersection(
        a->beg_x, a->beg_y,
        a->end_x, a->end_y,
        b->beg_x, b->beg_y,
        b->end_x, b->end_y,
        &ox, &oy
      );

      if (hit) {
        next->next = a;
        a->head = next->head;

        structure_test(b->head);
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


  ToolKind nearest_tool_kind;
  // float to_nearest_tool = 0;
  Hand *nearest_hand = 0;
  WoodDispenser *nearest_wood_dispenser = 0;
  {
    float to_nearest_hand = 0;
    float to_nearest_wood_dispenser = 0;
    float x = env.mouse_x;
    float y = env.mouse_y;
    find_nearest_hand(x, y, &nearest_hand, &to_nearest_hand);
    find_nearest_wood_dispenser(x, y, &nearest_wood_dispenser, &to_nearest_wood_dispenser);
    nearest_tool_kind = (to_nearest_wood_dispenser < to_nearest_hand) ? ToolKind_WoodDispenser : ToolKind_Hand;
    // to_nearest_tool = (to_nearest_wood_dispenser < to_nearest_hand) ? to_nearest_wood_dispenser : to_nearest_hand;
  }
  for (int i = 0; i < ARR_LEN(state.hands); i++) {
    Hand *hand = state.hands + i;
    if (hand->stage == HandStage_Uninit) continue;
    if (nearest_tool_kind == ToolKind_Hand && nearest_hand == hand) continue;
    if (state.hand_selected == hand) continue;

    HandOut hando = {0};
    hand_out(hand, state.tick, &hando);
    draw_hand(&hando, clr_hand);
  }
  for (int i = 0; i < ARR_LEN(state.wood_dispensers); i++) {
    WoodDispenser *wd = state.wood_dispensers + i;
    if (wd->stage == WoodDispenserStage_Uninit) continue;
    if (nearest_tool_kind == ToolKind_WoodDispenser && nearest_wood_dispenser == wd) continue;

    WoodDispenserOut wdo = {0};
    wood_dispenser_out(wd, state.tick, &wdo);
    draw_wood_dispenser(&wdo);
  }

  /* draw boards (all of 'em!) */
  for (int i = 0; i < ARR_LEN(state.boards); i++) {
    Board *b = state.boards + i;
    if (b->stage == BoardStage_Uninit) continue;

    if (b->stage == BoardStage_Fading) {
      float duration = b->tick_stage_end - b->tick_stage_start;
      float elapsed = state.tick - b->tick_stage_start;
      float t = (duration > 0) ? (elapsed / duration) : 0;
      rendr.alpha = ease_out_quad(1 - t);
    }

    plot_line_thick(b->beg_x, b->beg_y,
                    b->end_x, b->end_y,
                    clr_board, 0.18f);
    rendr.alpha = 1.0f;
  }

  /* {
    float x = env.mouse_x;
    float y = env.mouse_y;
    Color red = { 255, 0, 0, 255 };
    plot_line(x-0.1, y-0.0, x+0.1, y+0.0, red);
    plot_line(x-0.0, y-0.1, x+0.0, y+0.1, red);
  } */

  /* decorate house */
  for (int i = 0; i < ARR_LEN(state.houses); i++) {
    House *h = state.houses + i;
    if (h->stage == HouseStage_Uninit) continue;

    Board *b = h->boards + 0;
    Board *o = h->boards + 1;

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

    quad_to_tri(points, tri);

    if (h->stage == HouseStage_Init) {
      h->stage = HouseStage_Falling;
    } 
    if (h->stage == HouseStage_Falling) {
      Board *bs = h->boards + 0;

      /* center of base of triangle */
      float cx = (tri[1].x + tri[2].x) / 2;
      float cy = (tri[1].y + tri[2].y) / 2;
      float dx = tri[0].x - cx;
      float dy = tri[0].y - cy;

      float rot_now = atan2f(dy, dx);
      float ideal_rot = MATH_PI/2.0f;
      float delta_rot = rads_dist(rot_now, ideal_rot);
      board_pivot(bs+0, cx, cy, delta_rot*0.04f);
      board_pivot(bs+1, cx, cy, delta_rot*0.04f);

      bs[0].beg_y -= (cy - (cy*0.95f));
      bs[0].end_y -= (cy - (cy*0.95f));
      bs[1].beg_y -= (cy - (cy*0.95f));
      bs[1].end_y -= (cy - (cy*0.95f));

      if (delta_rot < 0.01f && cy <= 0.01f)
        h->stage = HouseStage_Settled;
    }

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

  ui(UiEvent_Render);

  state_save();
}

static void ui(UiEvent event) {
  int mousedown = event == UiEvent_Mousedown || event == UiEvent_DoubleClick;
  int double_click = event == UiEvent_DoubleClick;
  int render = event == UiEvent_Render;

  float x = env.mouse_x;
  float y = env.mouse_y;

  ToolKind nearest_tool_kind;
  float to_nearest_tool = 0;
  Hand *nearest_hand = 0;
  WoodDispenser *nearest_wood_dispenser = 0;
  {
    float to_nearest_hand = 0;
    float to_nearest_wood_dispenser = 0;
    find_nearest_hand(x, y, &nearest_hand, &to_nearest_hand);
    find_nearest_wood_dispenser(x, y, &nearest_wood_dispenser, &to_nearest_wood_dispenser);
    nearest_tool_kind = (to_nearest_wood_dispenser < to_nearest_hand) ? ToolKind_WoodDispenser : ToolKind_Hand;
    to_nearest_tool = (to_nearest_wood_dispenser < to_nearest_hand) ? to_nearest_wood_dispenser : to_nearest_hand;
  }

  int override = 0;
  switch (state.mode) {
    case (Mode_View): {
      if (to_nearest_tool < CLICK_DIST_SELECT) {
        if (nearest_tool_kind == ToolKind_Hand) {
          if (double_click) {
            state.hand_selected = nearest_hand;
            state.mode = Mode_HandView;
          }
          if (render) {
            override = 1;
            HandOut hando = {0};
            hand_out(nearest_hand, state.tick, &hando);
            draw_hand(&hando, (Color) { 185,  75,  55, 255 });
          }
        }
        if (nearest_tool_kind == ToolKind_WoodDispenser) {
          if (double_click) {
            state.wood_dispenser_selected = nearest_wood_dispenser;
            state.mode = Mode_WoodDispenserView;
          }
          if (render) {
            override = 1;
            WoodDispenserOut wdo = {0};
            wood_dispenser_out(nearest_wood_dispenser, state.tick, &wdo);
            rendr.alpha = 0.5f;
            draw_wood_dispenser(&wdo);
            rendr.alpha = 1.0f;
          }
        }
      }
    } break;

    case (Mode_HandView): {
      Color body_clr = { 215,  50,   0, 255 };

      Hand *hand = state.hand_selected;
      float to_selected_hand = mag(hand->x - env.mouse_x,
                                   hand->y - env.mouse_y);

      float    to_angle_grab = mag((hand->x + cosf(hand->   angle_grab)) - x,
                                   (hand->y + sinf(hand->   angle_grab)) - y);
      if (to_angle_grab < CLICK_DIST_SELECT) {
        if (mousedown) {
          state.mode = Mode_HandMoveAngleGrab;
          break;
        }
      }

      float to_angle_release = mag((hand->x + cosf(hand->      angle_release)) - x,
                                   (hand->y + sinf(hand->      angle_release)) - y);
      if (to_angle_release < CLICK_DIST_SELECT) {
        if (mousedown) {
          state.mode = Mode_HandMoveAngleRelease;
          break;
        }
      }

      if (to_selected_hand > CLICK_DIST_UNSELECT) {
        if (mousedown) {
          state.mode = Mode_View;
          state.hand_selected = 0;
        }
        if (render) 
          body_clr = (Color) { 180,  85,   0, 255 };
      }
      if (to_selected_hand < CLICK_DIST_SELECT) {
        if (mousedown) {
          state.mode = Mode_HandMove;
        }
        if (render)
          body_clr = (Color) { 255,   0,   0, 255 };
      }

      if (render) {
        HandOut hando = {0};
        hand_out(hand, state.tick, &hando);

        float x = hando.x;
        float y = hando.y;
        plot_line(x-0.1, y-0.0, x+0.1, y+0.0, body_clr);
        plot_line(x-0.0, y-0.1, x+0.0, y+0.1, body_clr);

        HandOut arc = hando;
        Color arc_clr = { 205, 140, 205, 255 };
        arc.theta = hando.angle_grab;
        arc.grip_width = GRIP_WIDTH_GRAB;
        draw_hand(&arc, arc_clr);
        {
          Color clr = arc_clr;
          float _x = x + cosf(arc.angle_grab);
          float _y = y + sinf(arc.angle_grab);
          if (mag(_x - env.mouse_x, _y - env.mouse_y) < CLICK_DIST_SELECT)
            clr.g = 0, clr.b = 0;
          plot_line(_x-0.1, _y-0.0, _x+0.1, _y+0.0, clr);
          plot_line(_x-0.0, _y-0.1, _x+0.0, _y+0.1, clr);
        }

        arc.theta = hando.angle_release;
        arc.grip_width = GRIP_WIDTH_RELEASE;
        draw_hand(&arc, arc_clr);
        {
          Color clr = arc_clr;
          float _x = x + cosf(arc.angle_release);
          float _y = y + sinf(arc.angle_release);
          if (mag(_x - env.mouse_x, _y - env.mouse_y) < CLICK_DIST_SELECT)
            clr.g = 0, clr.b = 0;
          plot_line(_x-0.1, _y-0.0, _x+0.1, _y+0.0, clr);
          plot_line(_x-0.0, _y-0.1, _x+0.0, _y+0.1, clr);
        }

        // override = 0;
        /* we want to override only to avoid us getting drawn twice */
        if (nearest_tool_kind == ToolKind_Hand &&
            nearest_hand == state.hand_selected )
          override = 1;
        draw_hand(&hando, body_clr);
      }

    } break;

    case (Mode_WoodDispenserView): {
      WoodDispenser *wood_dispenser = state.wood_dispenser_selected;
      float to_selected_wood_dispenser = mag(wood_dispenser->x - env.mouse_x,
                                             wood_dispenser->y - env.mouse_y);

      if (to_selected_wood_dispenser > CLICK_DIST_UNSELECT) {
        if (mousedown) {
          state.mode = Mode_View;
          state.wood_dispenser_selected = 0;
        }
      }
      if (to_selected_wood_dispenser < CLICK_DIST_SELECT) {
        if (mousedown) {
          state.mode = Mode_WoodDispenserMove;
        }
      }
      if (render) {
        float x = wood_dispenser->x;
        float y = wood_dispenser->y;
        plot_line(x-0.1, y-0.0, x+0.1, y+0.0, (Color) { 255, 0, 0, 255 });
        plot_line(x-0.0, y-0.1, x+0.0, y+0.1, (Color) { 255, 0, 0, 255 });
      }
    } break;

    case (Mode_BuyPreview): {
      if (state.preview_tool_kind == ToolKind_WoodDispenser) {
        if (render) {
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
        }
      }
      if (state.preview_tool_kind == ToolKind_Hand) {
        if (render) {
          HandOut ghost = {0};
          hand_out(
            &(Hand) { .x = env.mouse_x,
                      .y = env.mouse_y,
                      .tick_stage_start = state.tick-TICK_SECOND,
                      .tick_stage_end   = state.tick+1,
                      .stage = HandStage_Grab },
            state.tick,
            &ghost
          );
          draw_hand(&ghost, (Color) { 120, 185, 120, 255 });
        }
      }
    } break;

    case (Mode_HandMoveAngleGrab): {
      Hand *hand = state.hand_selected;
      HandOut hando = {0};
      hand_out(hand, state.tick, &hando);

      float x = hando.x;
      float y = hando.y;
      plot_line(x-0.1, y-0.0, x+0.1, y+0.0, (Color) { 255, 0, 0, 255 });
      plot_line(x-0.0, y-0.1, x+0.0, y+0.1, (Color) { 255, 0, 0, 255 });

      HandOut arc = hando;

      {
        Color arc_clr = { 140, 205, 185, 255 };
        arc.theta = atan2f(env.mouse_y - y, env.mouse_x - x);
        arc.grip_width = GRIP_WIDTH_GRAB;
        draw_hand(&arc, arc_clr);
      }

      {
        Color arc_clr = { 205, 140, 205, 255 };
        arc.theta = hando.angle_release;
        arc.grip_width = GRIP_WIDTH_RELEASE;
        draw_hand(&arc, arc_clr);
      }
    } break;

    case (Mode_HandMoveAngleRelease): {
      Hand *hand = state.hand_selected;
      HandOut hando = {0};
      hand_out(hand, state.tick, &hando);

      float x = hando.x;
      float y = hando.y;
      plot_line(x-0.1, y-0.0, x+0.1, y+0.0, (Color) { 255, 0, 0, 255 });
      plot_line(x-0.0, y-0.1, x+0.0, y+0.1, (Color) { 255, 0, 0, 255 });

      HandOut arc = hando;

      {
        Color arc_clr = { 205, 140, 205, 255 };
        arc.theta = hando.angle_grab;
        arc.grip_width = GRIP_WIDTH_GRAB;
        draw_hand(&arc, arc_clr);
      }

      {
        Color arc_clr = { 140, 205, 185, 255 };
        arc.theta = atan2f(env.mouse_y - y, env.mouse_x - x);
        arc.grip_width = GRIP_WIDTH_RELEASE;
        draw_hand(&arc, arc_clr);
      }
    } break;

    case (Mode_WoodDispenserMove): {
      WoodDispenserOut wdo = {0};
      wood_dispenser_out(nearest_wood_dispenser, state.tick, &wdo);
      rendr.alpha = 0.5;
      wdo.x = env.mouse_x;
      wdo.y = env.mouse_y;
      draw_wood_dispenser(&wdo);
      rendr.alpha = 1.0;
    } break;

    case (Mode_HandMove): {
      Hand *hand = state.hand_selected;
      HandOut hando = {0};
      hand_out(hand, state.tick, &hando);

      Color body_clr = {  80, 195,  80, 255 };

      HandOut arc = hando;
      arc.x = env.mouse_x;
      arc.y = env.mouse_y;

      override = 1;
      draw_hand(&hando, body_clr);
      draw_hand(&arc, body_clr);

      Color arc_clr = { 140, 205, 185, 255 };
      arc.theta = hando.angle_grab;
      arc.grip_width = GRIP_WIDTH_GRAB;
      draw_hand(&arc, arc_clr);

      arc.theta = hando.angle_release;
      arc.grip_width = GRIP_WIDTH_RELEASE;
      draw_hand(&arc, arc_clr);
    } break;
  }

  if (!override) {
    if (nearest_hand != 0 &&
        nearest_tool_kind == ToolKind_Hand) {
      HandOut hando = {0};
      hand_out(nearest_hand, state.tick, &hando);
      draw_hand(&hando, clr_hand);
    }
    if (nearest_wood_dispenser != 0 &&
        nearest_tool_kind == ToolKind_WoodDispenser) {
      WoodDispenserOut wdo = {0};
      wood_dispenser_out(nearest_wood_dispenser, state.tick, &wdo);
      draw_wood_dispenser(&wdo);
    }
  }
  /*
    switch (state.mode) {

      case (Mode_HandMove): {
        if (hand != state.hand_selected) break;
      } break;

      case (Mode_WoodDispenserView): {} break;
      case (Mode_WoodDispenserMove): {} break;
    }

    -- WOOD DISPENSER --

    if (state.mode == Mode_View)
      if (wd == near_mouse_wood_dispenser && to_near_mouse_wood_dispenser < CLICK_DIST_SELECT)
        rendr.alpha = 0.5f;
    if (state.mode == Mode_WoodDispenserView) {
      float x = wd->x;
      float y = wd->y;
      plot_line(x-0.1, y-0.0, x+0.1, y+0.0, (Color) { 255, 0, 0, 255 });
      plot_line(x-0.0, y-0.1, x+0.0, y+0.1, (Color) { 255, 0, 0, 255 });

      if (to_near_mouse_wood_dispenser < CLICK_DIST_SELECT)
        rendr.alpha = 0.30f;
      else
        rendr.alpha = 0.75f;
    }
    if (state.mode == Mode_WoodDispenserMove) {
      rendr.alpha = 0.35f;
      WoodDispenserOut wdo = {0};
      wood_dispenser_out(wd, state.tick, &wdo);
      wdo.x = env.mouse_x;
      wdo.y = env.mouse_y;
      draw_wood_dispenser(&wdo);
      rendr.alpha = 1.00f;
    } */

  {
    /* Math.ceil(Math.log10(Math.pow(2, 32))) */
    char buf[10] = {0};

    uint32_t n = state.tick;
    int nchars = 0;
    for (int i = n; i > 10; i /= 10) nchars++;

    for (int i = 0; i <= nchars; i++) {
      buf[nchars - i] = '0' + (n % 10);
      n /= 10;
    }
    TextPopup tp = {
      .msg = buf,
      .size = 0.058f,
      .y = 0.95f,
      .x = 0.01f
    };
    text_popup_draw(&tp);
  }

  {
    TextPopup tp = {
      .msg = "wood dispenser",
      .size = 0.028f,
      .y = 0.95f,
      .x = 0.78f
    };
    text_popup_draw(&tp);
    if (text_popup_hovered(&tp)) {
      text_popup_underline(&tp);

      if (mousedown)
        state.mode = Mode_BuyPreview,
        state.preview_tool_kind = ToolKind_WoodDispenser;
    }

    if (!text_popup_hovered(&tp)      &&
        mousedown &&
        state.mode == Mode_BuyPreview &&
        state.preview_tool_kind == ToolKind_WoodDispenser
    ) {

      WoodDispenser *wd = wood_dispenser_alloc();
      if (!wd) return;

      wd->stage = WoodDispenserStage_Init;
      wd->x = env.mouse_x;
      wd->y = env.mouse_y;

      state.mode = Mode_View;
    }
  }

  {
    TextPopup tp = {
      .msg = "hand",
      .size = 0.028f,
      .y = 0.95f,
      .x = 0.78f
    };
    tp.y -= tp.size * 2.2f;

    text_popup_draw(&tp);
    if (text_popup_hovered(&tp)) {
      text_popup_underline(&tp);

      if (mousedown) state.mode = Mode_BuyPreview,
                     state.preview_tool_kind = ToolKind_Hand;
    }

    if (!text_popup_hovered(&tp)      &&
        mousedown &&
        state.mode == Mode_BuyPreview &&
        state.preview_tool_kind == ToolKind_Hand
    ) {
      Hand *hand = hand_alloc();
      if (!hand) return;

      hand->stage = HandStage_Init;
      hand->x = env.mouse_x;
      hand->y = env.mouse_y;

      state.mode = Mode_View;
    }
  }
}
