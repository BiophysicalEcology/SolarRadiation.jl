import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}
const navTemp = {
  nav: [
    { text: 'Home', link: '/' },
    { text: 'Get Started', link: '/get_started' },
    { text: 'Manual',
      items: [
        { text: 'Introduction', link: '/manual/introduction' },
        { text: 'Solar geometry', link: '/manual/solar_geometry' },
        { text: 'Atmosphere', link: '/manual/atmosphere' },
        { text: 'Diffuse models', link: '/manual/diffuse_models' },
        { text: 'Terrain', link: '/manual/terrain' },
        { text: 'Poles, tropics and equator', link: '/manual/latitude_examples' },
        { text: 'Aerosols', link: '/manual/aerosols' },
        { text: 'Dates and calendars', link: '/manual/dates_calendars' },
        { text: 'Performance', link: '/manual/performance' },
        { text: 'Data tables', link: '/manual/data_tables' },
        { text: 'References', link: '/manual/references' },
      ]
    },
    { text: 'Tutorials',
      items: [
        { text: 'Mapping with GADS', link: '/tutorials/mapping' },
        { text: 'Elevation: the Himalaya', link: '/tutorials/elevation' },
        { text: 'Terrain: Saba', link: '/tutorials/saba' },
        { text: 'Photosynthetically active radiation', link: '/tutorials/par' },
        { text: 'Ultraviolet radiation and sunburn', link: '/tutorials/sunburn' },
      ]
    },
    { text: 'API', link: '/api' }
  ],
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: "Clear-sky spectral solar radiation for biophysical ecology",
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  ignoreDeadLinks: false,
  vite: {
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    build: {
      assetsInlineLimit: 0, // so we can tell whether we have created inlined images or not, we don't let vite inline them
    },
    optimizeDeps: {
      exclude: [
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ],
    },
    ssr: {
      noExternal: [
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ],
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },

  themeConfig: {
    outline: 'deep',
    // https://vitepress.dev/reference/default-theme-config
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
      { text: 'Get Started', link: '/get_started' },
      { text: 'Manual',
        items: [
          { text: 'Introduction', link: '/manual/introduction' },
          { text: 'Solar geometry', link: '/manual/solar_geometry' },
          { text: 'Atmosphere', link: '/manual/atmosphere' },
          { text: 'Diffuse models', link: '/manual/diffuse_models' },
          { text: 'Terrain', link: '/manual/terrain' },
          { text: 'Poles, tropics and equator', link: '/manual/latitude_examples' },
          { text: 'Aerosols', link: '/manual/aerosols' },
          { text: 'Dates and calendars', link: '/manual/dates_calendars' },
          { text: 'Performance', link: '/manual/performance' },
          { text: 'Data tables', link: '/manual/data_tables' },
          { text: 'References', link: '/manual/references' },
        ]
      },
      { text: 'Tutorials',
        items: [
          { text: 'Mapping with GADS', link: '/tutorials/mapping' },
          { text: 'Elevation: the Himalaya', link: '/tutorials/elevation' },
          { text: 'Terrain: Saba', link: '/tutorials/saba' },
          { text: 'Photosynthetically active radiation', link: '/tutorials/par' },
          { text: 'Ultraviolet radiation and sunburn', link: '/tutorials/sunburn' },
        ]
      },
      { text: 'API', link: '/api' }
    ],
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      // { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    }
  }
})
