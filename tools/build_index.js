const path = require('path')
const fs = require('fs')
const lunr = require('lunr')
const { parse } = require('node-html-parser')

const HTML_FOLDER = '../build/html'
const SEARCH_FIELDS = ['title', 'body']
const EXCLUDE_FILES = [
  '404.html',
  'search.html',
  'genindex.html',
  'tutorial-slides.html'
]
const OUTPUT_INDEX = 'lunr_index'

/**
 * Is HTML
 * @param {string} filename File name
 * @returns {boolean} HTML
 */
const isHtml = (filename) => {
  const lower = filename.toLowerCase()
  return lower.endsWith('.htm') || lower.endsWith('.html')
}

/**
 * Find HTML files
 * @param {string} folder Folder
 * @returns {Array} HTML files
 */
const findHtml = (folder) => {
  if (!fs.existsSync(folder)) {
    console.error('Could not find folder: ', folder)
    return []
  }

  const files = fs.readdirSync(folder, { withFileTypes: true })

  const htmls = []
  for (let file of files) {
    const fileName = file.name
    if (file.isDirectory()) {
      const dirName = path.join(folder, file.name)
      const otherHtmls = findHtml(dirName)
      htmls.push(...otherHtmls.map((o) => path.join(fileName, o)))
    } else {
      if (isHtml(fileName) && !EXCLUDE_FILES.includes(fileName))
        htmls.push(fileName)
    }
  }

  return htmls
}

/**
 * Read HTML
 * @param {string} root Root path
 * @param {string} file File name
 * @returns {Object} Data { link, title, body }
 */
const readHtml = (root, file) => {
  const fileName = path.join(root, file)
  const html = fs.readFileSync(fileName).toString()

  // Parse HTML
  const document = parse(html, {
    script: false,
    noscript: true,
    style: false,
    pre: true
  })

  // Title
  const titles = document.getElementsByTagName('title')
  const title = titles.find((t) => t.parentNode.rawTagName === 'head').text

  // Body
  const toctree = document.querySelector('.toctree-wrapper')
  document.removeChild(toctree)

  let body = document.querySelector('#content').text
  body = body.replaceAll('\\(', '').replaceAll('\\)', '').replaceAll('\\', '')
  body = body.replaceAll(/<\/?span[^>]*>/g, '')

  // Return
  return {
    link: file,
    title,
    body
  }
}

/**
 * Parse text to anchor
 * @param {string} text Text
 * @returns {string} Anchor
 */
const parseAnchor = (text) => {
  // First pass
  text = text.toLowerCase()
  text = text.trim()
  text = text.replace(/ /g, '-')
  text = text.replace(/[’&\/\\#,+()$~%.'":*?<>{}]/g, '-')

  // Second pass
  while (text.slice(-1) === '-') text = text.slice(0, -1)
  text = text.replace(/-+/g, '-')

  return text
}

/**
 * Read single HTML file
 * @param {string} root Root path
 * @param {string} file File name
 * @returns {Array} Data [{ id, link, title, body }, ...]
 */
const readSingleHtml = (root, file) => {
  const fileName = path.join(root, file)
  const html = fs.readFileSync(fileName).toString()

  // Parse HTML
  const document = parse(html, {
    script: false,
    noscript: true,
    style: false,
    pre: true
  })

  // Remove toctree
  const toctree = document.querySelector('.toctree-wrapper')
  document.removeChild(toctree)

  // Content
  const content = document.querySelector('#content')

  // Sections
  const sections = content.getElementsByTagName('section')

  return sections.map((section, index) => {
    let titles = section.getElementsByTagName('h1')
    if (!titles.length) titles = section.getElementsByTagName('h1')
    if (!titles.length) titles = section.getElementsByTagName('h2')
    if (!titles.length) titles = section.getElementsByTagName('h3')
    if (!titles.length) titles = section.getElementsByTagName('h4')
    if (!titles.length) titles = section.getElementsByTagName('h5')
    if (!titles.length) titles = section.getElementsByTagName('h6')

    const title = titles[0].text

    let body = section.text
    body = body.replaceAll('\\(', '').replaceAll('\\)', '').replaceAll('\\', '')
    body = body.replaceAll(/<\/?span[^>]*>/g, '')

    return {
      id: index,
      link: file + '#' + parseAnchor(title),
      title,
      body
    }
  })
}

/**
 * Build Lunr index
 * @param {Array} docs Docs
 * @returns {Array} Index
 */
const buildIndex = (docs) => {
  return lunr(function () {
    SEARCH_FIELDS.forEach((field) => {
      this.field(field)
    })
    docs.forEach((doc) => {
      this.add(doc)
    }, this)
  })
}

/**
 * Build preview
 * @param {Array} docs Docs
 * @returns {Array} Previews
 */
function buildPreviews(docs) {
  return docs.map((doc) => ({
    ...doc,
    preview: doc.body
  }))
}

/**
 * Main
 */
const main = () => {
  const files = findHtml(HTML_FOLDER)

  console.info('Building index for these files:')
  const docs = []
  const index_i = []
  const previews_i = []
  for (let file of files) {
    console.info('\t - ' + file)
    docs.push({
      id: docs.length,
      ...readHtml(HTML_FOLDER, file)
    })

    const doc = readSingleHtml(HTML_FOLDER, file)
    index_i.push(buildIndex(doc))
    previews_i.push(buildPreviews(doc))
  }
  const index = buildIndex(docs)
  const previews = buildPreviews(docs)

  let js = 'const LUNR_DATA = ' + JSON.stringify(index) + ';\n'

  js += 'const LUNR_PAGEDATA = ['
  index_i.forEach((idx_i, i) => {
    js += JSON.stringify(idx_i)
    if (i < index_i.length - 1) js += ',\n'
  })
  js += '];\n'

  js += 'const PREVIEW_DATA = ' + JSON.stringify(previews) + ';\n'

  js += 'const PREVIEW_PAGEDATA = ['
  previews_i.forEach((preview_i, i) => {
    js += JSON.stringify(preview_i)
    if (i < previews_i.length - 1) js += ',\n'
  })
  js += '];'

  fs.writeFile(OUTPUT_INDEX + '.js', js, (err) => {
    if (err) {
      return console.error(err)
    }
    console.info('Index saved as ' + OUTPUT_INDEX)
  })
}

main()
